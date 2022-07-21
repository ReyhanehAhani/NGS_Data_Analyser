import itertools
import shutil
import argparse
from multiprocessing.pool import ThreadPool as Pool

FILTER_INTRONIC = True

class ThreadedParser():
    def __init__(self, path): 
        self.path = path
        self.initColumns()

    def getNextLine(self):
        with open(self.path) as f:
            for line in f:
                yield line

    def initColumns(self): 
        self.columns = [x for x in next(self.getNextLine()).split(',')]

    def processLine(self, line): 
        pass

    def afterProcess(self): 
        pass


    def run(self):
        pool = Pool(processes=8)

        for line in self.getNextLine():
            pool.map(self.processLine, (line, ))
        
        pool.close()
        pool.join()

        return self.afterProcess(), self.columns


class ChildParser(ThreadedParser):
    childGenes = []
    sharedGenes = []

    def processLine(self, line):
        line = [x.strip('"') for x in line.split(',')]
        data = dict(zip(self.columns, line))
        
        if data['Hom Iranome'] != '.' and data['Hom Iranome'] != '0':
            return

        if data['Het Iranome'].isnumeric(): 
            if int(data['Het Iranome']) >= 80:
                return

        if data['Het Our DB'].isnumeric(): 
            if int(data['Het Our DB']) >= 40:
                return

        if FILTER_INTRONIC:
            if data['Func.refGene'] == 'intronic': 
                return
        
        if data['Zygosity'] == 'hom':
            return

        self.childGenes.append(data)

    def afterProcess(self):
        for key, group in itertools.groupby(self.childGenes, lambda x: x['Gene.refGene']):
            duplicate_group, group = itertools.tee(group)
            if len(tuple(duplicate_group)) > 1: 
                for item in group:
                    self.sharedGenes.append(item)

        return self.sharedGenes


class MotherParser(ThreadedParser):

    def __init__(self, path, childGenes):
        super().__init__(path)
        self.childGenes = childGenes
        self.sharedGenes = []

    motherChildSharedGenes = []

    def processLine(self, line):
        line = [x.strip('"') for x in line.split(',')]
        data = dict(zip(self.columns, line))

        if FILTER_INTRONIC:
            if data['Func.refGene'] == 'intronic': 
                return

        if data['Zygosity'] == 'hom':
            return
    
        for gene in self.childGenes:
            if gene['Start'] == data['Start'] and gene['End'] == data['End'] and gene['Gene.refGene'] == data['Gene.refGene']:
                self.sharedGenes.append(gene)

    def afterProcess(self):
        return self.sharedGenes

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--keep-intronic', action='store_true')
    parser.add_argument('--no-keep-intronic', dest='keep-intronic', action='store_false')
    parser.set_defaults(keep_intronic=False)
    parser.add_argument('mother', help="The file address for the mother gene csv file")
    parser.add_argument('child', help="The file address for the child gene csv file")

    args = parser.parse_args()

    childName = ''.join(args.child.split('.')[:-1])
    shutil.copyfile(args.child, f"{childName}.copy.csv")

    FILTER_INTRONIC = not args.keep_intronic

    childParser = ChildParser(args.child)
    childGenes, childColumns = childParser.run()  

    motherParser = MotherParser(args.mother, childGenes)
    sharedGenes = motherParser.run()

    with open(args.child, 'w') as f:
        f.write(','.join(childColumns) + '\n')
        
        for genes in sharedGenes:
            f.write(','.join(genes.values()) + '\n')
