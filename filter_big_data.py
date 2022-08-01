import itertools
import shutil
import argparse
from multiprocessing.pool import ThreadPool as Pool
import xlsxwriter
import time
import gzip

FILTER_INTRONIC = True

class ThreadedParser():
    def __init__(self, path): 
        self.path = path
        self.initColumns()

    def getNextLine(self):
        f = open(self.path)
        
        if self.path.endswith('.gz'):
            f = gzip.open(self.path, 'rt') 

        for line in f:
            yield line

    def initColumns(self): 
        self.columns = [x.strip('"') for x in next(self.getNextLine()).split(',')]

        indexer = lambda e: self.columns.index(e) if e in self.columns else None

        self.HOM_IRANOME = indexer('Hom Iranome')
        self.HET_IRANOME = indexer('Het Iranome')
        self.HET_OUR_DB = indexer('Het Our DB')
        self.FUNC_REFGENE = indexer('Func.refGene')
        self.ZYGOSITY = indexer('Zygosity')
        self.GENE_REFGENE = indexer('Gene.refGene')
        self.START = indexer('Start')
        self.END = indexer('End')


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
        data = [x.strip('"') for x in line.split(',')]
        
        if data[self.HOM_IRANOME] != '.' and data[self.HOM_IRANOME] != '0':
            return

        if data[self.HET_IRANOME].isnumeric(): 
            if int(data[self.HET_IRANOME]) >= 80:
                return

        if data[self.HET_OUR_DB].isnumeric(): 
            if int(data[self.HET_OUR_DB]) >= 40:
                return

        if FILTER_INTRONIC:
            if data[self.FUNC_REFGENE] == 'intronic': 
                return
        
        if data[self.ZYGOSITY] == "hom":
            return

        self.childGenes.append(data)

    def afterProcess(self):
        for key, group in itertools.groupby(self.childGenes, lambda x: x[self.GENE_REFGENE]):
            duplicate_group, group = itertools.tee(group)
            if len(tuple(duplicate_group)) > 1: 
                for item in group:
                    self.sharedGenes.append(item)
        
        return self.sharedGenes


class MotherFatherParser(ThreadedParser):
    def __init__(self, path, childGenes, childColumns):
        self.childGenes = childGenes
        self.childColumns = childColumns
        self.sharedGenes = []
        super().__init__(path)

    motherChildSharedGenes = []

    def initColumns(self):
        super().initColumns()
        
        child_indexer = lambda e: self.childColumns.index(e) if e in self.childColumns else None
        self.CHILD_GENE_REFGENE = child_indexer('Gene.refGene')
        self.CHILD_START = child_indexer('Start')
        self.CHILD_END = child_indexer('End')
        
    def processLine(self, line):
        data = [x.strip('"') for x in line.split(',')]

        if FILTER_INTRONIC:
            if data[self.FUNC_REFGENE] == 'intronic': 
                return

        if data[self.ZYGOSITY] == 'hom':
            return
    
        for gene in self.childGenes:
            if gene[self.CHILD_START] == data[self.START] and gene[self.CHILD_END] == data[self.END] and gene[self.CHILD_GENE_REFGENE] == data[self.GENE_REFGENE]:
                self.sharedGenes.append(gene)

    def afterProcess(self):
        return self.sharedGenes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--keep-intronic', action='store_true')
    parser.add_argument('--no-keep-intronic', dest='keep-intronic', action='store_false')
    parser.set_defaults(keep_intronic=False)

    subparsers = parser.add_subparsers(dest='mode', required=True)

    mother_child = subparsers.add_parser('mother_child', help='Filter the mother child dataset')
    mother_child.add_argument('mother', help="The file address for the mother gene csv file")
    mother_child.add_argument('child', help="The file address for the child gene csv file")

    father_mother_child = subparsers.add_parser('father_mother_child', help='Filter the father mother child dataset')
    father_mother_child.add_argument('father', help="The file address for the father gene csv file")
    father_mother_child.add_argument('mother', help="The file address for the mother gene csv file")
    father_mother_child.add_argument('child', help="The file address for the child gene csv file")
    
    args = parser.parse_args()

    childName = ''.join(args.child.split('.')[:-1])

    FILTER_INTRONIC = not args.keep_intronic

    childParser = ChildParser(args.child)
    childGenes, childColumns = childParser.run()
    print(f'Child filter done, found {len(childGenes)} genes')

    motherParser = MotherFatherParser(args.mother, childGenes, childColumns)
    sharedGenes, _ = motherParser.run()
    print(f'Mother filter done, found {len(sharedGenes)} genes')
    label = 'احتمال کامپوند در مادر و فرزند'

    if args.mode == 'father_mother_child':
        motherParser = MotherFatherParser(args.father, sharedGenes, childColumns)
        sharedGenes, _ = motherParser.run()
        print(f'Father filter done, found {len(sharedGenes)} genes')
        label = 'موارد کامپوند در فرزند'

    current_row = 0

    workbook = xlsxwriter.Workbook(f'{childName}.xlsx')
    worksheet = workbook.add_worksheet()
    merge_format = workbook.add_format({
        'bold':     True,
        'align':    'center',
        'valign':   'vcenter',
        'fg_color': 'black',
        'font_color': 'white',
        'font_size': 16
    })

    worksheet.write_row(0, 0, childColumns)
    current_row += 1

    worksheet.merge_range(current_row, 4, current_row, 7, label, merge_format)
    current_row += 1

    for genes in sharedGenes:
        worksheet.write_row(current_row, 0, genes)
        current_row += 1

    workbook.close()


if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f"--- {round(time.time() - start_time, 3)} seconds ---")
