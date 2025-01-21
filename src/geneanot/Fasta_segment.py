"""
Efficient Fasta segment retrival.
"""
import gzip as gz
import os
from typing import Callable, Generator


def myopen(file: str) -> Callable:
    """Support for both .gz open and text open."""
    return gz.open if str(file)[-3:] == '.gz' else open


class Fasta_segment():
    '''
    Efficient reads of a segments starting from a given offset
    from a FASTA file.

    This class can be used with Fasta files where ALL sequence
    rows (excluding headers) are of the same length!

    In case of a FASTA with multiple sequences (multiple headers), the header name
    of the corresponding sequence must be provided. The header name is the string
    that proceeds the character ">".
    For example, the header name of the header:

    ">NC_000011.10 Homo sapiens chromosome 11, GRCh38.p13 Primary Assembly"

    is "NC_000011.10".

    In case of a FASTA with a single sequence (single header), the header name
    is not required.

    The class reads only the segment into memory and not the all FASTA sequence.
    This can be used, for example, to read a segment from a FASTA file
    that contains a chromosome of the human genome.
    '''
    def __init__(self, files_info=None) -> None:
        '''
        files_info (optional input) is a dictionary:
        For a FASTA with multiple headers, the key is a fasta file name and
        the header name as follows: "<fasta file name>:<header name>".
        For a FASTA with a single header the key is simply the file name (the header
        can also be provided, in which case the key will include it as well, however this is not required).
        The value is a dictionary that contains basic statistics of the fasta file.
        Keys of files_info[file] are:
            'init_loc': offset of the first NT in file (0 indicated the first NT)
            'l_size': number of NTs per line excluding \n,
            'line_size': number of NTs per line including \n.
            'file_size': total number of NTs in file, excluding \n (and excluding the header). This is
                         available only for FASTA files with a single header.

            For example: files_info = {'file1.fas': {'init_loc': 70, 'l_size': 80, 'line_size': 81, 'file_size': 240},
                               'file2.fas:NC_000011.10': {'init_loc': 70, 'l_size': 80, 'line_size': 81}}

        '''
        self.files_info = files_info if files_info else {}

    def __repr__(self) -> str:
        return "Fasta_segment([files_info])"

    def __str__(self) -> str:
        return "Efficient reads of segments from FASTA files.\n"

    def __len__(self) -> int:
        return len(self.files_info)

    def show_processed_files(self) -> list:
        return list(self.files_info.keys())

    def show_files_info(self) -> None:
        print(self.str_files_info())

    def str_files_info(self) -> str:
        s = ''
        for k, v in self.files_info.items():
            s += f"{k}:\n"
            for ki, vi in v.items():
                s += f"\t{ki}: {vi}\n"

        return s

    def read_segment(self, file: str, offset: int, size: int) -> str:
        '''
        Reads a segment from a FASTA with a single sequence (single header).
        file - FASTA file name
        offset - 0 based offset of the start segment in file (from the first sequence NT). 
                 For example, value of 1 indicates starting from the second NT in the file.
        size - number of NTs in a segment to read.
        '''
        if not os.path.exists(file):
            print(f"Error: file {file} does not exists !!")
            return ""

        with myopen(file)(file, 'rt') as fp:
            if file in self.files_info:
                l_size = self.files_info[file]['l_size']
                line_size = self.files_info[file]['line_size']
                init_loc = self.files_info[file]['init_loc']
            else:
                init_loc, line_size, l_size = self.__compute_file_stats(fp)
                self.files_info[file] = {'init_loc': init_loc, 'l_size': l_size, 'line_size': line_size}

            #offset_num_lines, offset_in_line = offset // l_size, np.mod(offset, l_size)
            offset_num_lines, offset_in_line = divmod(offset, l_size)
            # accounting for extra characters due to multiple \n that will be discarded
            add_for_newline = ((offset+size-1) // l_size) - offset_num_lines
            fp.seek(init_loc + offset_num_lines*line_size + offset_in_line)
            return fp.read(size+add_for_newline).replace('\n', '')

    def __compute_file_stats(self, fp) -> tuple:
        '''
        Computs the stats that are needed to read a segment from
        a FASTA file with a single sequence.
        '''
        # account for a header (if exists) in the first line of the file
        fp.seek(0, os.SEEK_SET)
        if fp.read(1) == ">":
            fp.readline()
            init_loc = fp.tell()
        else:
            init_loc = 0

        # init_loc is the location (0-based) of the first NT in file
        fp.seek(init_loc)
        fp.readline()  # advance handle to begining of second line
        line_size = fp.tell() - init_loc # number of NTs per line, including the \n
        return init_loc, line_size, line_size-1

    def total_num_chars(self, file: str) -> int:
        '''
        Returns the total number of characters (NTs) in a FASTA 
        file (excluding the header and \n) without reading the all file.

        This function supports only Fasta with a single headr.
        Support for FASTAs with multiple headers is TBD.
        '''
        if file in self.files_info:
            return self.files_info[file]['file_size']
        return self.__compute_total_num_NTs(file)

    def __compute_total_num_NTs(self, file: str) -> int:
        '''
        Computes the total number of characters (NTs) in a FASTA
        file (excluding the header and \n) and updates self.files_info.
        '''
        with myopen(file)(file, 'rt') as fp:
        #with open(file, 'rt') as fp:
            init_loc, line_size, l_size = self.__compute_file_stats(fp)

            # check if last character is new line
            last_loc = fp.seek(0, os.SEEK_END)
            fp.seek(fp.tell() - 1, os.SEEK_SET) 
            last_char = fp.read()

        delta = last_loc-init_loc
        q, r = divmod(delta, line_size)
        num_chars = delta - q # remove counts of new lines
        if q==1 and r==0 and (last_char!='\n'):
            num_chars += 1 # if only one line in file and last character is not newline add one

        self.files_info[file] = {'init_loc': init_loc,
                                 'l_size': l_size,
                                 'line_size': line_size,
                                 'file_size': num_chars}

        return num_chars

    def multiple_headers_read_segment(self, file: str, header_name: str, offset: int, size: int) -> str:
        '''
        Reads from a FASTA file with multiple sequences (multiple headers).
        The offset and size here are relative to the sequence with header name header_name.

        offset - a 0-based offset from the beginning of the sequence following a header with a header name
                 header_name.

        For example, for header_name='name1', offset=3, size=4, the return sequence is of size 4 NTs, starting
        from (the 0-based) offset 3 of the sequence that follows the header with header name header_name.

        Recall that, for example, the header name of the header:

        ">NC_000011.10 Homo sapiens chromosome 11, GRCh38.p13 Primary Assembly"

        is "NC_000011.10".
        '''
        init_loc = None
        token = file + ':' + header_name
        with myopen(file)(file, 'rt') as fp:
            if token in self.files_info:
                l_size = self.files_info[token]['l_size']
                line_size = self.files_info[token]['line_size']
                init_loc = self.files_info[token]['init_loc']
            else:
                line = fp.readline()
                while line:
                    if line[0] == '>':
                        # header name
                        name = line.split()[0][1:]
                        if name == header_name:
                            init_loc = fp.tell()  # file offset of the first NT of the corresponding sequence
                            break
                    line = fp.readline()
                if init_loc is None:
                    print(f"Did not find header name {header_name} in file {file} !!")
                    return ''
                fp.seek(init_loc)
                fp.readline()
                line_size = fp.tell() - init_loc
                l_size = line_size - 1
                self.files_info[token] = {'init_loc': init_loc, 'l_size': l_size, 'line_size': line_size}

            #offset_num_lines, offset_in_line = offset // l_size, np.mod(offset, l_size)
            offset_num_lines, offset_in_line = divmod(offset, l_size)
            # accounting for extra characters due to multiple \n that will be discarded
            add_for_newline = ((offset+size-1) // l_size) - offset_num_lines
            fp.seek(init_loc + offset_num_lines*line_size + offset_in_line)
            return fp.read(size+add_for_newline).replace('\n', '')

    def fasta_gen(self, file: str, segment_size: int, init_offset: int=0, jump_size: int=0) -> Generator[str, None, None]:
        '''
        Fasta generator.
        Returns an iterator for reading segments of size segment_size NTs, separated by
        jump_size NTs, starting from an offset (0-based) init_offset.

        To read segments consecutively, set jump_size to 0 (which is the default value).
        '''
        try:
            with myopen(file)(file, 'rt') as fp:
            #with open(file, 'rt') as fp:
                # handeling the header
                fp.seek(0, os.SEEK_SET)
                if fp.read(1) == ">":
                    fp.readline()
                    init_loc = fp.tell()
                else:
                    init_loc = 0

                # computing number of characters (including \n) per line
                fp.seek(init_loc)
                fp.readline()
                line_size = fp.tell() - init_loc # number of NTs per line, including the \n
                ext_newlines = segment_size // line_size # number of \n in a segment_size read with offset 0 
                jump_newlines = jump_size // line_size  # same for jump_size

                # starting from init_offset (accounting for \n)
                fp.seek(init_loc+init_offset+(init_offset//line_size))

                while True:
                    segment = fp.read(segment_size+ext_newlines).replace('\n', '')
                    # might need another single read as ext_newlines assumed reading from beginning of the line
                    if len(segment)<segment_size:
                        segment += fp.read(1).replace('\n', '')

                    if segment != '':
                        yield segment.replace('\n', '')
                    else:
                        break #raise StopIteration

                    # jump to next segment (the jump is excluding \n)
                    if fp.read(jump_size+jump_newlines).count('\n') > jump_newlines:
                        fp.seek(fp.tell()+1)

        except IOError as err:
            print(f"Error: {file} not found !! ({err})")
