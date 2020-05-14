'''
Utility classes
'''
import logging
import os
import subprocess


class SysUtil(object):

    @classmethod
    def mkdir_p(cls, dir_path):
        assert isinstance(dir_path, str), dir_path
        if dir_path == '':
            return
        if os.path.exists(dir_path) and os.path.isdir(dir_path):
            return
        try:
            os.makedirs(dir_path)
        except Exception:
            f'makedirs({dir_path}) failed'
            raise
        assert os.path.isdir(dir_path), f'{dir_path} is not a directory'

    @classmethod
    def is_executable_file(cls, file_path):
        if not os.path.isfile(file_path):
            return False
        if not os.access(file_path, os.X_OK):
            return False
        return True

    @classmethod
    def run_cmd(cls, cmd_line, timeout=None, dry_run=False):
        """
        Run `cmd_line` in shell with 'timeout'.
        If `dry_run` is True, only print the command line to screen
        """
        assert isinstance(cmd_line, str)
        assert isinstance(dry_run, bool)
        logging.info(f'running command line: "{cmd_line}" with timeout={timeout}, dry_run={dry_run}')
        if dry_run:
            print(cmd_line)
            return
        try:
            completed_proc = subprocess.run(cmd_line, shell=True, check=True, timeout=timeout,
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.SubprocessError as err:
            logging.info(f'command line {cmd_line} with timeout={timeout} failed with error:\n{err}')
            print(f'command line {cmd_line} with timeout={timeout} failed with error\n{err}')
            raise
        logging.info(f'command line run finished, output will be returned')
        std_out = completed_proc.stdout.decode()
        return std_out
