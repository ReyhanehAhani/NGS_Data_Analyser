import argparse
from multiprocessing.pool import ThreadPool as Pool
import gzip
import pathlib
import matplotlib.pyplot as plt

from typing import Union, Any, List

import base64

Path = Union[str, pathlib.Path]

class ThreadedParser():
    def __init__(self, path: Path): 
        self.path = path

    def getNextLine(self):
        f = open(self.path)
        
        if self.path.endswith('.gz'):
            f = gzip.open(self.path, 'rt') 

        for line in f:
            yield line

    def processLine(self, line: str): 
        pass

    def afterProcess(self) -> Any: 
        pass

    def run(self):
        pool = Pool(processes=8)

        for line in self.getNextLine():
            pool.map(self.processLine, (line, ))
        
        pool.close()
        pool.join()

        return self.afterProcess()

class MetaDataParser(ThreadedParser):
    def __init__(self, path: Path, phenotypes: List[str]):
        super().__init__(path)
        self.found = dict.fromkeys(phenotypes)
        for k in self.found.keys():
            self.found[k] = []

    def processLine(self, line: str):
        data = tuple(x.strip() for x in line.split('\t'))

        for phenotype in self.found.keys():
            if phenotype.lower() in line.lower():
                self.found[phenotype].append(data[0])
                
    def afterProcess(self):
        return self.found

class ChrParser(ThreadedParser):
    def __init__(self, path: Path, chr: str, pos: str, ref: str, alt: str):
        super().__init__(path)

        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt

        self.found = {}

    def processLine(self, line):
        data = tuple(x.strip() for x in line.split('\t'))

        if data[0] != self.chr:
            return

        if data[1] != self.pos:
            return

        if data[2] != self.ref:
            return

        found = False    
        for alt in data[3].split(','):
            if alt == self.alt:
                found = True
                break

        if not found:
            return

        alts = tuple(x for x in data[3].split(',') if x)
        zygosities = tuple(x for x in data[4].split(';') if x)
        ids = tuple(x for x in data[5].split(';') if x)
        
        # 0/1: Het for first alt
        # 1/1: Hom for first alt

        # 0/2: Het for second alt
        # 2/2: Hom for second alt

        # 0/3: Het for third alt
        # 3/3: Hom for third alt

        for zygosity, id in zip(zygosities, ids): 
            self.found[id] = 'het' if '0' in zygosity else 'hom'

    def afterProcess(self): 
        return self.found

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata', help="Path to metadata file")
    parser.add_argument('chromosomes', help="Path to chromosomes file")
    parser.add_argument('chr', help="Snip chromosome type")
    parser.add_argument('pos', help="Snip position")
    parser.add_argument('ref', help="Snip reference type")
    parser.add_argument('alt', help="Snip alternative type")
    parser.add_argument('phenotypes', help="Pheneotypes to search")
    args = parser.parse_args()

    phenotypes = MetaDataParser('TestData_metaData.txt', args.phenotypes.split(',')).run()
    print('Found Phenotypes')
    chromosomes = ChrParser('chr11Data_test.txt.gz', args.chr, args.pos, args.ref, args.alt).run()
    print('Found chromosomes')

    fig, ax = plt.subplots()
    plt.grid(zorder=0)
    
    max_freq = 0

    for index, (phenotype, pids) in enumerate(phenotypes.items()): 
        nHet = 0
        nHom = 0

        for pid in pids: 
            if pid in chromosomes.keys():
                if chromosomes[pid] == 'het':
                    nHet += 1
                else:
                    nHom += 1
        
        if len(pids) > 0:
            freq = (1 * nHet + 2 * nHom) / (2 * len(pids))
            
            if freq > max_freq:
                max_freq = freq
        else: 
            freq = 0

        text = f'''Freq: {freq}\nn(Het): {nHet}\nn(Hom): {nHom}'''
        ax.bar(index, freq, width=1, edgecolor="white", linewidth=0.7)
        ax.text(index - 0.25, 0.05 + freq, text)
        
    ax.set_ylim((0, max_freq + 0.25))
    ax.set_xlabel('Phenotypes')
    ax.set_ylabel('Frequency')
    ax.set_xticks(range(len(phenotypes.items())), labels=phenotypes.keys(), rotation=45)

    figure = plt.gcf()
    #figure.subplots_adjust(bottom=0.4)
    figure.set_size_inches(20, 18)
    plt.savefig(f"{args.chr}_{args.pos}_{args.ref}_{args.alt}.png", dpi=100)

if __name__ == "__main__":
    main()
