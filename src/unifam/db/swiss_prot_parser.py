'''
SwissProt .dat parser
'''
import re


class ID(object):
    def __init__(self, entry_name, status, length):
        assert isinstance(entry_name, str), entry_name
        assert isinstance(status, str), status
        assert isinstance(length, int), length
        self.entry_name = entry_name
        self.status = status
        self.length = length


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
        self._name_list = name_list
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


class SwissProtParser(object):
    '''
    Class to parse SwissProt .dat (annotation) flat file.
    We will only implement the annotation fields that unifam cares about.
    '''

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
        full_name_list = cls._get_value_for_key('Full')
        short_name_list = cls._get_value_for_key('Short')
        ec_list = cls._get_value_for_key('EC')
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
        assert line.startswith('GN '), line
        name_list = cls._get_value_for_key('Name')
        synonym_list = cls._get_value_for_key('Synonyms')
        ORF_list = cls._get_value_for_key('ORFNames')
        OLN_list = cls._get_value_for_key('OrderedLocusNames')
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

        format:
        OX   Taxonomy_database_Qualifier=Taxonomic code;

        Example
        -------
        OX   NCBI_TaxID=9606;
        OX   NCBI_TaxID=562;

        Returns
        -------
        '''
        assert line.startswith('OX '), line
        tax_ids = cls._get_value_for_key(line, 'NCBI_TaxID')
        assert len(tax_ids) == 1, f'length of tax_ids is not 1: {tax_ids}'
        return tax_ids[0]
