'''
SwissProt .dat parser
'''
import re
import os
import json
from unifam.base.util import SysUtil
from collections import defaultdict
from collections import ChainMap
from collections import Counter


class ClusterAnnot(object):
    '''
    Class to store annotation for a protein cluster that contain at least 1 swiss prot protein.

    Parameters
    ----------
    annot_dict_list : list of dict
        Each dict corresponds to the annotation of a protein in the cluster.
        Its keys should be a super set of ProteinAnnot.CLUSTER_ANNOT_KEYS.
    '''

    def __init__(self, annot_dict_list):
        assert isinstance(annot_dict_list, list)
        assert annot_dict_list
        for annot_dict in annot_dict_list:
            assert isinstance(annot_dict, dict)
            assert set(annot_dict) >= set(ProteinAnnot.CLUSTER_ANNOT_KEYS), \
                    f'annot_dict missing keys: {set(ProteinAnnot.CLUSTER_ANNOT_KEYS) - set(annot_dict)}'
        self._annot_dict_list = annot_dict_list

    @classmethod
    def from_prot_annot_list(cls, prot_annot_list):
        '''
        class method to create instance from list of ProteinAnnot instances
        '''
        assert isinstance(prot_annot_list, list), type(prot_annot_list)
        assert prot_annot_list
        prot_id_list = [prot_annot.get_id() for prot_annot in prot_annot_list]
        assert len(prot_id_list) == len(set(prot_id_list)), f'prot_annot_list has duplicate proteins'
        annot_dict_list = [prot_annot.to_cluster_annot_dict() for prot_annot in prot_annot_list]
        return cls(annot_dict_list)

    def __len__(self):
        return len(self._annot_dict_list)

    def get_cluster_annot_dict(self, share_thresh=None):
        '''
        When there is more than 1 protein in the cluster,
        an annotation value should only be used for the cluster
        if it's shared by >= `share_thresh` * len(cluster) proteins.

        Parameters
        ----------
        share_thresh : real
            real number in [0,1]
        '''
        if share_thresh is None:
            share_thresh = 0.5
        assert share_thresh >= 0. and share_thresh <= 1., share_thresh

        annot_dict = dict()
        for k in ProteinAnnot.CLUSTER_ANNOT_KEYS:
            annot_dict[k] = self._agg_annot(k, share_thresh)
        return annot_dict

    def _agg_annot(self, annot_key, share_thresh):
        # concatenate all the annot value lists together to a big list
        annot_val_list = sum([prot_annot_dict[annot_key]
                              for prot_annot_dict in self._annot_dict_list],
                             [])
        val_counter = Counter(annot_val_list)
        count_thresh = share_thresh * len(self)
        return [annot_val for annot_val, count in val_counter.items()
                if count >= count_thresh]


class ProteinAnnot(object):
    '''
    class to store annotation for a single protein
    '''
    CLUSTER_ANNOT_KEYS = ['kw', 'GO', 'name', 'ORF', 'OLN', 'synonym', 'full_name', 'short_name', 'ec']

    def __init__(self, prot_id, accession_number_list, prot_description,
                 prot_gene_name, organism, lineage,
                 ncbi_tax_id, kw_list, cross_ref_dict):
        assert isinstance(prot_id, ID)
        assert isinstance(prot_description, DE)
        assert isinstance(prot_gene_name, GN)
        assert isinstance(accession_number_list, list)
        assert isinstance(organism, str)
        assert isinstance(lineage, list)
        assert isinstance(ncbi_tax_id, str)
        assert isinstance(kw_list, list)
        assert isinstance(cross_ref_dict, dict)
        self.prot_id = prot_id
        self.accession_number_list = accession_number_list
        self.prot_description = prot_description
        self.prot_gene_name = prot_gene_name
        self.organism = organism
        self.lineage = lineage
        self.ncbi_tax_id = ncbi_tax_id
        self.kw_list = kw_list
        self.cross_ref_dict = cross_ref_dict
        self._annot_dict = dict(ChainMap(
            self.prot_id.to_dict(),
            self.prot_description.to_dict(),
            self.prot_gene_name.to_dict(),
            {'accession': self.accession_number_list,
             'organism': self.organism,
             'lineage': self.lineage,
             'ncbi_tax_id': self.ncbi_tax_id,
             'kw': self.kw_list,
             'GO': self.cross_ref_dict.get('GO', []),
            }))

    def __eq__(self, other):
        assert isinstance(other, self.__class__)
        return self.prot_id.entry_name == other.prot_id.entry_name

    def get_id(self):
        return self.prot_id.entry_name

    def to_dict(self):
        return self._annot_dict.copy()

    def get_pathologic_dict(self):
        '''
        Returns a dictionary with annotations that will be used in pathway tools.
        Some fields are not from annotation, such as STARTBASE and ENDBASE.

        Fields in Pathologic (.pf) file for a gene:
        ABUNDANCE optional
        CODING-SEGMENT optional
        DBLINK optional <DB:Accession> (Use UNIPROT or SP for UniProt accession number)
        EC recommended
        STARTBASE, ENDBASE optional
        FUNCTION required, use "ORF" if unknown
        FUNCTION-SYNONYM optional
        GENE-COMMENT optional
        GO recommended (will be written in DBLINK field)
        ID highly recommended - unique identifier for the gene
        METACYC optional - MetaCyc reaction ID assigned for the protein if known
        NAME required - common name of the gene
        PRODUCT-TYPE required - P, PSEUDO, TRNA, RRNA, MISC-RNA
        SYNONYM

        Example
        -------
        ;;; The PF file format.
        ;;; This is a comment in front of an imaginary example record.
        ;;; It starts with a semicolon. Each record starts with an "ID"
        ;;; line, and is terminated by a "//" line.
        ID b1262
        NAME trpC
        STARTBASE 1317812
        ENDBASE 1316451
        DBLINK SP:P00909
        PRODUCT-TYPE P
        SYNONYM foo
        SYNONYM foo2
        GENE-COMMENT f453; 99 pct identical to TRPC_ECOLI SW:P00909
        ;;; The following shows how information about multiple functions of
        ;;; a protein is supplied:
        FUNCTION N-(5-phosphoribosyl) anthranilate isomerase
        EC 5.3.1.24
        FUNCTION-SYNONYM phosphoribosyl anthranilate isomerase
        FUNCTION-COMMENT Amino acid biosynthesis: Tryptophan (3rd step)
        FUNCTION indole-3-glycerolphosphate synthetase
        EC 4.1.1.48
        FUNCTION-COMMENT	Amino acid biosynthesis: Tryptophan (4th step)
        DBLINK GO:0000250
        //
        '''
        patho_dict = dict()
        db_link_list = [f'SP:{ac_num}' for ac_num in self.accession_number_list]
        if 'GO' in self.cross_ref_dict:
            db_link_list += self.cross_ref_dict['GO']
        patho_dict['DBLINK'] = db_link_list

        patho_dict['EC'] = self.prot_description.get_ec_list()
        patho_dict['FUNCTION'] = (self.prot_description.get_full_name_list() +
                                  self.prot_description.get_short_name_list())
        if not len(patho_dict['FUNCTION']):
            patho_dict['FUNCTION'] = ['ORF']
        patho_dict['FUNCTION-SYNONYM'] = self.prot_gene_name.get_synonym_list()
        patho_dict['FUNCTION-COMMENT'] = (self.prot_gene_name.get_ORF_list() +
                                          self.prot_gene_name.get_OLN_list())
        patho_dict['NAME'] = self.prot_gene_name.get_name_list()
        patho_dict['PRODUCT-TYPE'] = ['P']

        return patho_dict

    def to_cluster_annot_dict(self):
        """
        Returns a dict containing a subset of annotation information,
        that can be used as a protein cluster's annotation.
        """
        return {key: val for key, val in self._annot_dict.items() if key in self.CLUSTER_ANNOT_KEYS}


class ID(object):
    def __init__(self, entry_name, status, length):
        assert isinstance(entry_name, str), entry_name
        assert isinstance(status, str), status
        assert isinstance(length, int), length
        self.entry_name = entry_name
        self.status = status
        self.length = length

    def to_dict(self):
        return {'entry_name': self.entry_name,
                'status': self.status,
                'length': self.length}


class DE(object):
    def __init__(self, full_name_list=None, short_name_list=None, ec_list=None):
        if full_name_list is None:
            full_name_list = []
        if short_name_list is None:
            short_name_list = []
        if ec_list is None:
            ec_list = []
        assert isinstance(full_name_list, list)
        assert isinstance(short_name_list, list)
        assert isinstance(ec_list, list)
        self._full_name_list = full_name_list
        self._short_name_list = short_name_list
        self._ec_list = ec_list

    def get_full_name_list(self):
        return self._full_name_list.copy()

    def get_short_name_list(self):
        return self._short_name_list.copy()

    def get_ec_list(self):
        return self._ec_list.copy()

    def __add__(self, other_de):
        assert isinstance(other_de, self.__class__)
        return self.__class__(self._full_name_list + other_de._full_name_list,
                              self._short_name_list + other_de._short_name_list,
                              self._ec_list + other_de._ec_list)

    def to_dict(self):
        return {'full_name': self.get_full_name_list(),
                'short_name': self.get_short_name_list(),
                'ec': self.get_ec_list()}


class GN(object):
    '''
    OLN_list : list of str
        OrderedLocusNames
    '''

    def __init__(self, name_list=None, synonym_list=None, ORF_list=None,
                OLN_list=None):
        if name_list is None:
            name_list = []
        if synonym_list is None:
            synonym_list = []
        if ORF_list is None:
            ORF_list = []
        if OLN_list is None:
            OLN_list = []
        assert isinstance(name_list, list)
        assert isinstance(synonym_list, list)
        assert isinstance(ORF_list, list)
        assert isinstance(OLN_list, list)
        self._name_list = [name.upper() for name in name_list]
        self._synonym_list = synonym_list
        self._ORF_list = ORF_list
        self._OLN_list = OLN_list

    def get_name_list(self):
        return self._name_list.copy()

    def get_synonym_list(self):
        return self._synonym_list.copy()

    def get_ORF_list(self):
        return self._ORF_list.copy()

    def get_OLN_list(self):
        return self._OLN_list.copy()

    def __add__(self, other_de):
        assert isinstance(other_de, self.__class__)
        return self.__class__(self._name_list + other_de._name_list,
                              self._synonym_list + other_de._synonym_list,
                              self._ORF_list + other_de._ORF_list,
                              self._OLN_list + other_de._OLN_list)

    def to_dict(self):
        return {'name': self.get_name_list(),
                'ORF': self.get_ORF_list(),
                'OLN': self.get_OLN_list(),
                'synonym': self.get_synonym_list()}


class SwissProtParser(object):
    '''
    Class to parse SwissProt .dat (annotation) flat file.
    We will only implement the annotation fields that unifam cares about.
    '''
    @classmethod
    def save_annot(cls, id_to_annot, save_path):
        '''
        Save the parsed protein annotation `id_to_annot` to a file at `save_path`.
        '''
        assert isinstance(id_to_annot, dict)
        SysUtil.mkdir_p(os.path.dirname(save_path))
        with open(save_path, 'w') as save_f:
            json.dump({prot_id: prot_annot.to_dict() for prot_id, prot_annot in id_to_annot.items()},
                      save_f, indent=2)

    @classmethod
    def load_annot(cls, load_path):
        '''
        Load the json dumped (by `save_annot`) file to protein annotation dict.
        '''
        assert os.path.isfile(load_path), load_path
        with open(load_path, 'r') as f:
            return json.load(f)

    @classmethod
    def read_annot_file(cls, annot_file, start_pos=0, max_records=-1):
        '''
        Read SwissProt annotation (flat) file to ProteinAnnot objects.

        Parameters
        ----------
        annot_file: str
            path to the SwissProt .dat file
        start_pos: int
            starting pos (in bytes) for the parsing
        max_records: int
            maximum number of records to parse before return
            Default value -1 means there is no limit.
        '''
        assert os.path.isfile(annot_file), f'{annot_file} does not exist'
        assert isinstance(start_pos, int) and start_pos >= 0, start_pos
        assert isinstance(max_records, int), max_records
        assert max_records == -1 or max_records > 0, max_records

        id_to_annot = dict()
        should_parse = True if start_pos == 0 else False
        with open(annot_file, 'r') as annot_f:
            annot_f.seek(start_pos)
            record_lines = []
            # Use readline so that we can use annot_f.tell()
            # When using loop `for line in annot_f`, annot_f.tell() cannot be used
            line = annot_f.readline()
            while line:
                if line[:2] == '//':
                    if should_parse:
                        protein_annot = cls.get_protein_annot(cls.record_lines_to_dict(record_lines))
                        id_to_annot[protein_annot.get_id()] = protein_annot
                        num_annot = len(id_to_annot)
                        if max_records != -1 and num_annot >= max_records:
                            print(f'Reached max: {num_annot} records, last read position: {annot_f.tell()}')
                            break
                        if num_annot % 1000 == 0:
                            print(f'Read {num_annot} records...')
                    else:
                        print(f'Found // at pos {annot_f.tell()}')
                        should_parse = True
                    record_lines = []
                else:
                    record_lines.append(line)
                line = annot_f.readline()
            last_pos = annot_f.tell()
        print(f'Read {num_annot} records, last read position: {last_pos}')
        return id_to_annot

    @classmethod
    def record_lines_to_dict(cls, record_lines):
        '''
        Convert list of record lines for a protein annotation
        to a dictionary: code --> str
        with only the interested codes as keys.

        Among the code_list, GN, KW, and DR codes are optional,
        other codes have at least 1 line of content.
        '''
        assert isinstance(record_lines, list), type(record_lines)
        code_list = ['ID', 'AC', 'DE', 'GN', 'OS', 'OC', 'OX', 'KW', 'DR']
        code_to_content = defaultdict(str)
        for line in record_lines:
            if line[0] == ' ':
                continue
            code, content = line.split(maxsplit=1)
            if code in code_list:
                if code_to_content[code]:
                    code_to_content[code] += content if code == 'DR' else content.strip('\n')
                else:
                    code_to_content[code] = line if code == 'DR' else line.strip('\n')
        return code_to_content

    @classmethod
    def get_protein_annot(cls, code_to_content):
        return ProteinAnnot(
            prot_id=cls.read_id_line(code_to_content['ID']),
            accession_number_list=cls.read_ac_line(code_to_content['AC']),
            prot_description=cls.read_de_line(code_to_content['DE']),
            prot_gene_name=cls.read_gn_line(code_to_content['GN']),
            organism=cls.read_os_line(code_to_content['OS']),
            lineage=cls.read_oc_line(code_to_content['OC']),
            ncbi_tax_id=cls.read_ox_line(code_to_content['OX']),
            kw_list=cls.read_kw_line(code_to_content['KW']),
            cross_ref_dict=cls.read_dr_line(code_to_content['DR'])
        )

    @classmethod
    def read_id_line(cls, line):
        '''
        IDentification line, 1 and only 1 line

        Example
        -------
        ID   EntryName Status; SequenceLength.
        '''
        assert line.startswith('ID '), line
        assert line.endswith(' AA.'), line
        cols = line.split()
        entry_name = cols[1]
        status = cols[2].strip(' ;')
        length = int(cols[3])
        return ID(entry_name, status, length)

    @classmethod
    def read_ac_line(cls, line):
        '''
        ACcession number line(s), at least 1 line

        Example
        -------
        AC   AC_number_1;[ AC_number_2;]...[ AC_number_N;]
        '''
        assert line.startswith('AC '), line
        ac_numbers = [ac_number.strip(';') for ac_number in line.split()[1:]]
        return ac_numbers

    @classmethod
    def _get_value_for_key(cls, line, key):
        '''
        From a line in sprot.dat, find the values for a given key
        Return list.
        Format: key=val1, val2 {evidence}, ..., valn {evidence};
        '''
        values = []
        start_pos = 0
        field_pos = line.find(f'{key}=', start_pos)
        while(field_pos != -1):
            end_pos = line.find(';', field_pos)
            value = line[(field_pos + len(key) + 1):end_pos]
            values.extend([x.strip() for x in cls._strip_evidence(value).split(',')])
            start_pos = end_pos + 1
            field_pos = line.find(f'{key}=', start_pos)
        return values

    @classmethod
    def _strip_evidence(cls, string):
        '''
        Strip the {evidence} string from the value string

        Example: val2 {evidence} -> val2
        '''
        return re.sub(r'{.*}', '', string.strip())

    @classmethod
    def read_de_line(cls, line):
        '''
        We only care about Full name, Short name, and EC numbers.
        DEscription line(s), at least should include 1 line with RecName

        Example
        -------
        DE   RecName: Full=Annexin A5;
        DE            Short=Annexin-5;
        DE   AltName: Full=Annexin V;
        ...
        '''
        assert line.startswith('DE '), line
        full_name_list = cls._get_value_for_key(line, 'Full')
        short_name_list = cls._get_value_for_key(line, 'Short')
        ec_list = cls._get_value_for_key(line, 'EC')
        return DE(full_name_list, short_name_list, ec_list)

    @classmethod
    def read_gn_line(cls, line):
        '''
        Gene Name line, optional

        Example
        -------
        GN   Name=<name>; Synonyms=<name1>[, <name2>...]; OrderedLocusNames=<name1>[, <name2>...];
        GN   ORFNames=<name1>[, <name2>...];
        '''
        if line == '':
            return GN()
        assert line.startswith('GN '), line
        name_list = cls._get_value_for_key(line, 'Name')
        synonym_list = cls._get_value_for_key(line, 'Synonyms')
        ORF_list = cls._get_value_for_key(line, 'ORFNames')
        OLN_list = cls._get_value_for_key(line, 'OrderedLocusNames')
        return GN(name_list, synonym_list, ORF_list, OLN_list)

    @classmethod
    def read_os_line(cls, line):
        '''
        OS (Organism Species) line specifies the organism
        which was the source of the stored sequence. Once or more.

        We only keep the latin name, i.e. the part before the first parenthesis.
        Also need to concatenate the lines that continue from previous lines:
            line should end with  "."

        Example
        -------
        OS   Avian leukosis virus RSA (RSV-SRA) (Rous sarcoma virus (strain
        OS   Schmidt-Ruppin A)).
        '''
        assert line.startswith('OS '), line
        assert line.endswith('.'), f'line not complete: {line}'
        organism_latin_name = line.split('(', maxsplit=1)[0][5:]
        return organism_latin_name

    @classmethod
    def read_oc_line(cls, line):
        '''
        OC (Organism Classification) lines contain the taxonomic classification of the source organism.
        Once or more.

        Example
        -------
        OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia;
        OC   Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae;
        OC   Homo.

        Returns
        -------
        list of str
            taxonomy sequence of the source organism
        '''
        assert line.startswith('OC '), line
        assert line.endswith('.'), f'line not complete: {line}'
        return [x.strip() for x in line[:-1].split(maxsplit=1)[-1].split(';')]

    @classmethod
    def read_ox_line(cls, line):
        '''
        The OX (Organism taxonomy cross-reference) line is used to indicate the identifier
        of a specific organism in a taxonomic database.
        Once.
        We only care about NCBI_TaxID here.

        format:
        OX   Taxonomy_database_Qualifier=Taxonomic code;

        Example
        -------
        OX   NCBI_TaxID=9606;
        OX   NCBI_TaxID=562;

        Returns
        -------
        ncbi_tax_id: str
        '''
        assert line.startswith('OX '), line
        tax_ids = cls._get_value_for_key(line, 'NCBI_TaxID')
        assert len(tax_ids) == 1, f'length of tax_ids is not 1: {tax_ids}'
        return tax_ids[0]

    @classmethod
    def read_kw_line(cls, line):
        '''
        The KW (KeyWord) lines provide information that can be used to generate
        indexes of the sequence entries based on functional, structural, or other categories.
        Optional

        format:
        KW   Keyword[; Keyword...].

        Example
        -------
        KW   Protease inhibitor; Proteoglycan; Serine protease inhibitor; Signal;
        KW   Transmembrane; Zinc.

        Returns
        -------
        list of str
            list of keywords
        '''
        if line == '':
            return []
        assert line.startswith('KW '), line
        return [x.strip() for x in line[:-1].split(maxsplit=1)[-1].split(';')]

    @classmethod
    def read_dr_line(cls, line):
        '''
        The DR (Database cross-Reference) lines are used as pointers to information in external data resources that is
        related to UniProtKB entries.
        Optional

        format:
        DR   RESOURCE_ABBREVIATION; RESOURCE_IDENTIFIER; OPTIONAL_INFORMATION_1\
                [; OPTIONAL_INFORMATION_2][;OPTIONAL_INFORMATION_3].

        Example
        -------
        DR   EMBL; AY548484; AAT09660.1; -; Genomic_DNA.
        DR   RefSeq; YP_031579.1; NC_005946.1.
        DR   SwissPalm; Q6GZX4; -.
        DR   GeneID; 2947773; -.
        DR   KEGG; vg:2947773; -.
        DR   Proteomes; UP000008770; Genome.
        DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.

        Returns
        -------
        dict: str --> list of str
            resource_abbr --> resource identifier
        '''
        if line == '':
            return dict()
        assert line.startswith('DR '), line
        cross_ref_dict = defaultdict(list)
        for resource_line in line.split(maxsplit=1)[-1].split('.\n'):
            field_list = [x.strip() for x in resource_line.split(';') if x.strip() != '-']
            if len(field_list) > 1:
                cross_ref_dict[field_list[0]].append(field_list[1])
        return dict(cross_ref_dict)
