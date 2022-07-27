from Bio.Emboss import PrimerSearch
import re
from typing import Sequence
from Bio.Seq import Seq
from fasta_reader import read_fasta, write_fasta


## We nee a way to represent 
class PrimerRecord():

    def name():
        doc = """name"""
        def fget(self):
            return self._name

        def fset(self, value):
            self._name = value

        def fdel(self):
            del self._name
        return locals()
    name = property(**name())


    def start_position():
        doc = """start position"""
        def fget(self):
            return self._start_position

        def fset(self, value):
            self._start_position = value

        def fdel(self):
            del self._start_position
        return locals()
    start_position = property(**start_position())

    def end_position():
        doc = """end position"""
        def fget(self):
            return self._end_position

        def fset(self, value):
            self._end_position = value

        def fdel(self):
            del self._end_position
        return locals()
    end_position = property(**end_position())

    def __str__(self):
        return f'{self.name}: {self.start_position} - {self.end_position}'


start_matcher = re.compile(r"at ([0-9]+)")
end_matcher = re.compile(r"at \[([0-9]+)\]")

    


def load_amplifiers():
    primer_records = {}
    search_results = None    
    with open('./somefile.primersearch') as f:
        search_results = PrimerSearch.read(f)

    # We're only using one primer sequence
    amplifiers = search_results.amplifiers
    amplifier = amplifiers[list(amplifiers.keys())[0]]

    
    for result in amplifier:
        lines = result.hit_info.split('\n')
        record = PrimerRecord()
        name = lines[0].strip()
        record.name = name
        start_pos = int(start_matcher.search(lines[2]).group(1))
        record.start_position = start_pos
        end_pos = int(end_matcher.search(lines[3]).group(1))
        record.end_position = end_pos
        primer_records[record.name] = record


    
    return primer_records



def generate_fasta_file(primer_records):
    output = write_fasta('lactic_acid_pcr.fa.gz')
    for item in read_fasta('./lactic_acid_only.fa.gz'):
        id = item.defline.split(' ')[0].strip()
        if id in primer_records:
            primer = primer_records[id]
            seq = item.sequence[primer.start_position:(primer.end_position + 1)]
            if len(seq) > 0:
                output.write_item(item.defline,seq)
        

primer_info = load_amplifiers()
generate_fasta_file(primer_info)








    
