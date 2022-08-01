import pandas as pd
import numpy as np 
import os
import argparse
import io
import pathlib

from typing import Callable, Tuple, List, Union, Dict

Path = Union[str, pathlib.Path]

def generate_file_name(*files: List[str]) -> str: 
    files = [file.split(os.sep)[-1].replace('filtered_', '').replace('.csv','')[:-2] for file in files]
    return f"filtered_{'&'.join(files)}.xlsx"

class GeneralParser():
    match_columns = ["Het Iranome", "Hom Iranome", "Het Our DB", "Chr", "Start", "End", "Ref", "Alt", "Zygosity", "Gene.refGene", "ExonicFunc.refGene"]
    gene_exceptions = ['frameshift insertion', 'frameshift deletion', 'stopgain', 'stoploss', 'splicing']
    
    def warn(self, message: str):
        print("Warning:", message)

    def error(self, message: str):
        print("Error:", message)

    def success(self, message: str):
        print(message)

    def read_OMIMfile(self, path: Path) -> Dict[str, str]: 
            raw_genes: Dict[str, List[str]] = {}
            filtered_genes: Dict[str, str] = {}

            with open(path) as f: 
                for line in f:
                    line = line.strip().split('\t')

                    if raw_genes.get(line[3]): 
                        raw_genes[line[3]].append(line[2].strip())
                    else: 
                        raw_genes[line[3]] = [line[2]]
            
            AD = ['AD', 'XL', 'XLD', 'Smu', 'Mu', 'SMo', 'IC', 'YL']
            AR = ['AR', 'XLR', 'DR']

            for key, items in raw_genes.items():
                for item in items: 
                    if any(x in item for x in AR): 
                        filtered_genes[key] = 'AR'
                        break
                    elif any(x in item for x in AD):
                        filtered_genes[key] = 'AD'
                        break

            return filtered_genes

    def save_xlsx(self, dataframes: List[Tuple[pd.DataFrame, str]], output: Path):
        current_row: int = 0
        writer = pd.ExcelWriter(output, engine='xlsxwriter')

        # Get the worksheet and create the style
        workbook = writer.book
        worksheet = workbook.add_worksheet('Sheet1')
        writer.sheets['Sheet1'] = worksheet

        merge_format = workbook.add_format({
            'bold':     True,
            'align':    'center',
            'valign':   'vcenter',
            'fg_color': 'black',
            'font_color': 'white',
            'font_size': 16
        })

        for index, (df, label) in enumerate(dataframes):
            if index == 0:
                columns = pd.DataFrame(columns=df.columns)
                columns.to_excel(writer, sheet_name='Sheet1', header=True, index=False, startrow=current_row)
                current_row += 1

            worksheet.merge_range(current_row, 4, current_row, 7, label, merge_format)
            current_row += 1
            df.sort_values('Gene.refGene').to_excel(writer, sheet_name='Sheet1', header=False, index=False, startrow=current_row)
            current_row += df.shape[0]

        writer.save()

    def read_faulty_csv(self, path: Path) -> pd.DataFrame: 
        df: pd.DataFrame | None = None
        columns: List[str] | None = None

        with open(path, 'r') as f:
            columns = list(dict.fromkeys([x for x in f.readline().strip().split(',') if x != '']))
            df = pd.DataFrame()

            data_buffer = []
            
            for index, line in enumerate(f): 
                print(index, end='\r')

                line = line.strip()
                data = line.split(',')[:len(columns)]
                data = [x.replace('"', '').strip() for x in data]
                data_buffer.append(data)

                if index % 250 == 0:
                    small_df = pd.DataFrame(data_buffer)
                    df = pd.concat([small_df, df], ignore_index=True)
                    data_buffer.clear()

            if data_buffer:
                small_df = pd.DataFrame(data_buffer)
                df = pd.concat([small_df, df], ignore_index=True)
                data_buffer.clear()

            df.columns = columns
        
        return df

    def read_csv(self, path: Path, parent: str, data_filter: Callable[[pd.DataFrame, str, bool], pd.DataFrame], keep_intronic: bool = False) -> pd.DataFrame:
        try:
            df = pd.read_csv(path, index_col=False)
        except:
            self.warn(f'{path} contains bad lines, trying slower method')
            df = self.read_faulty_csv(path)
        
        self.success(f'{path} was read successfully')

        return data_filter(df, parent, keep_intronic)

    def concat_dataframes(self, dataframes: List[pd.DataFrame]) -> pd.DataFrame:
        dataframes = [df.astype(str).reset_index(drop=True) for df in dataframes]
        df = pd.concat(dataframes, ignore_index=True).reset_index(drop=True).sort_values(['Gene.refGene'])

        df['_rank'] = df.groupby(['Gene.refGene', 'Parent']).cumcount()
        df.sort_values(['Gene.refGene', '_rank'], inplace=True)
        df.drop(labels=['_rank'], axis=1, inplace=True)

        return df

    def mother_and_father_share_gene(self, row: pd.Series) -> bool:
        row = row.sort_values('Parent')
        parents = tuple(row['Parent'].values)
        return ('mother' in parents) and ('father' in parents)
    
    def compound_gene(self, row: pd.Series, match_columns: List[str]) -> bool:
        row = row.sort_values('Parent')
        parents = tuple(row['Parent'].values)

        if len(parents) != 2:
            return False

        if not (('mother' in parents) and ('father' in parents)):
            return False

        father = row[row['Parent'] == 'father']
        mother = row[row['Parent'] == 'mother']
        
        return not np.array_equal(father[match_columns].values, mother[match_columns].values)

    def mother_and_child_share_gene(self, row: pd.Series, check_zygosity: bool = True) -> bool:
        row = row.sort_values('Parent')

        if check_zygosity:
            return (tuple(row['Parent'].values) == ('child', 'mother')) and (tuple(row['Zygosity'].values) == ('hom', 'het'))
        else: 
            return (tuple(row['Parent'].values) == ('child', 'mother'))

    def mother_father_and_child_do_not_share_gene(self, row: pd.Series)  -> bool:
        row = row.sort_values('Parent')
        return ('father' not in row['Parent'].values) and ('mother' not in row['Parent'].values)

    def not_shared_path(self, row: pd.Series): 
        row = row[row['Parent'] != 'child']
        return len(row['Parent'].values) == 1

    def mother_and_child_or_father_and_child_share_gene(self, row : pd.Series) -> bool:
        row = row.sort_values('Parent')
        zygosity = tuple(row['Zygosity'].values)
        parents = tuple(row['Parent'].values)
        if ('child' in parents and (('father' in parents) or ('mother' in parents))):
            if zygosity[0] == 'hom' and zygosity[1] == 'het':
                if len(zygosity) > 2:
                    if zygosity[2] == 'het':
                        return True
                    else:
                        False
                return True
        return False

    def mother_and_child_share_path(self, row: pd.Series, omim: Dict[str, str], omim_check: str = 'AR') -> bool: 
        row = row.sort_values(['Parent'])

        if tuple(row['Parent'].values) != ('child', 'mother'):
            return False

        if tuple(row['Zygosity'].values) != ('het', 'het'):
            return False
        
        gene = row['Gene.refGene'].values[0]

        if omim.get(gene, 'AR') == omim_check:
            return True
        else:
            return False

    def for_check_in_mother_child(self, row: pd.Series, omim: Dict[str, str], omim_check: str = 'AR'): 
        if row['Zygosity'] != 'hom':
            return False
        
        if omim.get(row['Gene.refGene'], 'AR') == omim_check:
            return True

        return False

    def for_check_in_father_mother_child(self, row: pd.Series, omim: Dict[str, str]): 
        if row['Zygosity'] == 'het' and omim.get(row['Gene.refGene'], 'AR') == 'AD':
            return True
        if row['Zygosity'] == 'hom' and omim.get(row['Gene.refGene'], 'AR') == 'AR':
            return True

        return False

    def filter_normal(self, df: pd.DataFrame, parent: str, keep_intronic: bool = False) -> pd.DataFrame:
        df = df[(df['Hom Iranome'] == '0') | (df['Hom Iranome'] == '.')]
        
        if not keep_intronic:
            df = df[df['Func.refGene'] != 'intronic']

        df['Het Iranome'] = pd.to_numeric(df['Het Iranome'], errors='coerce')
        df = df[(df['Het Iranome'].isna()) | (df['Het Iranome'] < 80)]
        df['Het Iranome'] = df['Het Iranome'].fillna('.')

        df['Het Our DB'] = pd.to_numeric(df['Het Our DB'], errors='coerce')
        df = df[(df['Het Our DB'].isna()) | (df['Het Our DB'] < 40)]
        df['Het Our DB'] = df['Het Our DB'].fillna('.')

        zygosity_index = df.columns.get_loc('Zygosity')
        df.insert(zygosity_index, 'Parent', parent)

        return df

    def filter_path(self, df: pd.DataFrame, parent: str, keep_intronic: bool = True) -> pd.DataFrame:

        df = df[(df['Hom Iranome'] == '0') | (df['Hom Iranome'] == '.')]
        df = df[~(df['CLNSIG'].str.contains('Benign', case=False))]

        df = df[((df['ValueInfo2'].str.split(':').str[0] == '0/1')
                | (df['ValueInfo2'].str.split(':').str[0] == '1/1'))
                & (df['ValueInfo2'].str.split(':').str[1].str.split(',').str[1].astype(int) > 7)]

        zygosity_index = df.columns.get_loc('Zygosity')
        df.insert(zygosity_index, 'Parent', parent)

        return df

    def run(self):
        raise NotImplementedError()


class FatherMotherParser(GeneralParser):
    def __init__(self, mother: Path, father: Path, 
                       mother_path: Path, father_path: Path, 
                       output: Path,
                       keep_intronic: bool = False):
        self.mother = mother
        self.father = father
        self.mother_path = mother_path
        self.father_path = father_path
        self.output = output
        self.keep_intronic = keep_intronic
    
    def run(self): 
        mother = self.read_csv(self.mother, 'mother', self.filter_normal, self.keep_intronic)
        father = self.read_csv(self.father, 'father', self.filter_normal, self.keep_intronic)

        mother_path = self.read_csv(self.mother_path, 'mother', self.filter_path)
        father_path = self.read_csv(self.father_path, 'father', self.filter_path)

        normal_df = self.concat_dataframes([mother, father])
        path_df = self.concat_dataframes([mother_path, father_path])

        mother_and_father_shared_gene = normal_df.groupby(self.match_columns).filter(self.mother_and_father_share_gene)
        normal_df.drop(labels=mother_and_father_shared_gene.index, inplace=True)

        compound_gene = normal_df.groupby('Gene.refGene').filter(lambda row: self.compound_gene(row, self.match_columns))
        normal_df.drop(labels=compound_gene.index, inplace=True)
        
        dangerous_gene = normal_df[(normal_df['ExonicFunc.refGene'].isin(self.gene_exceptions))
                          | (normal_df['ExonicFunc.ensGene'].isin(self.gene_exceptions))
                          | (normal_df['ExonicFunc.knownGene'].isin(self.gene_exceptions))
                          | (normal_df['Func.refGene'].isin(self.gene_exceptions))
                          | (normal_df['Function_description'].isin(self.gene_exceptions))]
        dangerous_gene = dangerous_gene[dangerous_gene['Func.refGene'] != 'intronic']
        normal_df.drop(labels=dangerous_gene.index, inplace=True)
        
        mother_and_father_shared_path = path_df.groupby(self.match_columns).filter(self.mother_and_father_share_gene)
        path_df.drop(labels=mother_and_father_shared_path.index, inplace=True)
        
        mother_and_father_not_shared_path = path_df[~(path_df.duplicated(subset=self.match_columns, keep=False))]
        path_df.drop(labels=mother_and_father_not_shared_path.index, inplace=True)
        
        mother_hom_check = path_df[(path_df['Zygosity'] == 'hom') & (path_df['Parent'] == 'mother')]
        path_df.drop(labels=mother_hom_check.index, inplace=True)
        
        father_hom_check = path_df[(path_df['Zygosity'] == 'hom') & (path_df['Parent'] == 'father')]
        path_df.drop(labels=father_hom_check.index, inplace=True)
        
        shared_genes = self.concat_dataframes([mother_and_father_shared_gene, mother_and_father_shared_path])
        
        datasets = [(shared_genes, 'موارد مشترک در زوج'),
                    (compound_gene, 'ژن مشترک برای احتمال کامپوند'),
                    (dangerous_gene, 'موارد خطرناک در هر یک از زوجین'),
                    (mother_and_father_not_shared_path, 'موارد پاتوژن غیرمشترک'),
                    (father_hom_check, 'برای بررسی در پدر'),
                    (mother_hom_check, 'برای بررسی در مادر')]

        self.save_xlsx(datasets, self.output)

class MotherChildParser(GeneralParser):
    def __init__(self, mother: Path, child: Path, 
                       mother_path: Path, child_path: Path, 
                       omim: Path,
                       output: Path,
                       keep_intronic: bool = False):
        self.mother = mother
        self.child = child
        self.mother_path = mother_path
        self.child_path = child_path
        self.output = output
        self.keep_intronic = keep_intronic
        self.omim = omim
    
    def run(self): 
        mother = self.read_csv(self.mother, 'mother', self.filter_normal, self.keep_intronic)
        child = self.read_csv(self.child, 'child', self.filter_normal, self.keep_intronic)

        mother_path = self.read_csv(self.mother_path, 'mother', self.filter_path)
        child_path = self.read_csv(self.child_path, 'child', self.filter_path)

        omim_file = self.read_OMIMfile(self.omim)

        normal_df = self.concat_dataframes([mother, child])
        path_df = self.concat_dataframes([mother_path, child_path])

        for_check = path_df.groupby(self.match_columns).filter(lambda x: self.mother_and_child_share_path(x, omim_file, 'AD'))        
        carrier_chance = path_df.groupby(self.match_columns).filter(lambda x: self.mother_and_child_share_path(x, omim_file, 'AR'))

        shared_mother_child_gene = normal_df.groupby(self.match_columns).filter(self.mother_and_child_share_gene)
        normal_df.drop(labels=shared_mother_child_gene.index, inplace=True)
        
        shared_mother_child_path = path_df.groupby(self.match_columns).filter(self.mother_and_child_share_gene)
        path_df.drop(labels=shared_mother_child_path.index, inplace=True)
        
        not_shared_mother_child_gene = normal_df.groupby('Gene.refGene').filter(self.mother_father_and_child_do_not_share_gene)
        normal_df.drop(labels=not_shared_mother_child_gene.index, inplace=True)
        
        not_shared_mother_child_path = path_df.groupby('Gene.refGene').filter(self.mother_father_and_child_do_not_share_gene)
        path_df.drop(labels=not_shared_mother_child_path.index, inplace=True)
        
        shared_mother_child_path_without_het = path_df.groupby(self.match_columns).filter(lambda x: self.mother_and_child_share_path(x, omim_file, 'AR'))
        path_df.drop(labels=shared_mother_child_path_without_het.index, inplace=True)
        

        shared_gene = self.concat_dataframes([shared_mother_child_gene, shared_mother_child_path])
        not_shared_path = pd.merge(path_df, self.concat_dataframes([shared_mother_child_path_without_het, shared_mother_child_path]), how='left', indicator=True).query("_merge == 'left_only'").drop('_merge', axis=1)

        dataframes = [
            (shared_gene, 'موارد مشترک در مادر و فرزند'),
            (not_shared_mother_child_gene, 'موارد غیرمشترک در فرزند (‌فقط فرزند)'),
            (not_shared_mother_child_path, 'موارد پاتوژن فرزند'),
            (mother_path[mother_path['Zygosity'] == 'hom'], 'موارد پاتوژن مادر'),
            (carrier_chance, 'احتمال ناقل بودن در مادر و فرزند'),
            (shared_mother_child_path_without_het, 'موارد مشترک پاتوژن مشترک در مادر و فرزند'),
            (for_check[for_check['Parent'] == 'child'], 'برای بررسی در فرزند'),
            (for_check[for_check['Parent'] == 'mother'], 'برای بررسی در مادر'),
            (not_shared_path, 'موارد غیرمشترک پاتوژن در مادر و فرزند')
        ]

        self.save_xlsx(dataframes, self.output)



class FatherMotherChildParser(GeneralParser): 
    def __init__(self, 
                 mother: Path, father: Path, child: Path,
                 mother_path: Path, father_path: Path, child_path: Path,
                 omim: Path,
                 output: Path,
                 keep_intronic: bool = False): 
        self.mother = mother
        self.father = father
        self.child = child
        self.mother_path = mother_path
        self.father_path = father_path
        self.child_path = child_path
        self.output = output
        self.keep_intronic = keep_intronic
        self.omim = omim

    def run(self): 
        father = self.read_csv(self.father, 'father', self.filter_normal, self.keep_intronic)
        mother = self.read_csv(self.mother, 'mother', self.filter_normal, self.keep_intronic)
        child = self.read_csv(self.child, 'child', self.filter_normal, self.keep_intronic)
        omim_file = self.read_OMIMfile(self.omim)

        father_path = self.read_csv(self.father_path, 'father', self.filter_path)
        mother_path = self.read_csv(self.mother_path, 'mother', self.filter_path)
        child_path = self.read_csv(self.child_path, 'child', self.filter_path)

        normal_df = self.concat_dataframes([father, mother, child])
        path_df = self.concat_dataframes([father_path, mother_path, child_path])

        for_check = path_df.groupby(self.match_columns).filter(lambda x: self.mother_and_child_share_path(x, omim_file, 'AD'))

        shared_gene = normal_df.groupby([x for x in self.match_columns if x != 'Zygosity']).filter(self.mother_and_child_or_father_and_child_share_gene)
        normal_df.drop(labels=shared_gene.index, inplace=True)
        
        shared_gene_path = path_df.groupby([x for x in self.match_columns if x != 'Zygosity']).filter(self.mother_and_child_or_father_and_child_share_gene)
        path_df.drop(labels=shared_gene_path.index, inplace=True)
        
        not_shared_gene = normal_df.groupby('Gene.refGene').filter(self.mother_father_and_child_do_not_share_gene)
        normal_df.drop(labels=not_shared_gene.index, inplace=True)
        
        not_shared_gene_path = path_df.groupby('Gene.refGene').filter(self.mother_father_and_child_do_not_share_gene)
        path_df.drop(labels=not_shared_gene_path.index, inplace=True)
        
        shared_gene_without_child = normal_df.groupby([x for x in self.match_columns if x != 'Zygosity']).filter(self.mother_and_father_share_gene)
        normal_df.drop(labels=shared_gene_without_child.index, inplace=True)
        
        shared_gene_without_child_path = path_df.groupby([x for x in self.match_columns if x != 'Zygosity']).filter(self.mother_and_father_share_gene)
        path_df.drop(labels=shared_gene_without_child_path.index, inplace=True)
        
        shared_gene_without_child = shared_gene_without_child[shared_gene_without_child['Parent'] != 'child']
        shared_gene_without_child_path = shared_gene_without_child_path[shared_gene_without_child_path['Parent'] != 'child']

        compound_gene = normal_df.groupby('Gene.refGene').filter(lambda row: self.compound_gene(row, self.match_columns))
        normal_df.drop(labels=compound_gene.index, inplace=True)
        
        father_mother_child_shared = self.concat_dataframes([shared_gene, shared_gene_path])
        father_mother_child_not_shared = self.concat_dataframes([not_shared_gene, not_shared_gene_path])
        father_mother_shared = self.concat_dataframes([shared_gene_without_child, shared_gene_without_child_path])

        dangerous_gene = normal_df[(normal_df['ExonicFunc.refGene'].isin(self.gene_exceptions))
                          | (normal_df['ExonicFunc.ensGene'].isin(self.gene_exceptions))
                          | (normal_df['ExonicFunc.knownGene'].isin(self.gene_exceptions))
                          | (normal_df['Func.refGene'].isin(self.gene_exceptions))
                          | (normal_df['Function_description'].isin(self.gene_exceptions))]
        
        dangerous_gene = dangerous_gene[dangerous_gene['Parent'] != 'child']
        dangerous_gene = dangerous_gene[dangerous_gene['Func.refGene'] != 'intronic']
        normal_df.drop(labels=dangerous_gene.index, inplace=True)
        
        mother_and_father_not_shared_path = path_df.groupby('Gene.refGene').filter(self.not_shared_path)
        mother_and_father_not_shared_path = mother_and_father_not_shared_path[mother_and_father_not_shared_path['Parent'] != 'child']
        path_df.drop(labels=mother_and_father_not_shared_path.index, inplace=True)

        
        not_shared_path = path_df

        datasets = [
            (father_mother_child_shared, 'موارد مشترک در پدر و مادر و فرزند'),
            (father_mother_child_not_shared, 'موارد غیرمشترک در پدر و مادر و فرزند'),
            (father_mother_shared, 'موارد مشترک در زوج'),
            (compound_gene, 'ژن مشترک برای احتمال کامپوند'),
            (dangerous_gene, 'موارد خطرناک در هر یک از زوجین'),
            (mother_and_father_not_shared_path, 'موارد پاتوژن غیر مشترک در زوج'),
            (for_check[for_check['Parent'] == 'father'], 'برای بررسی در پدر'),
            (for_check[for_check['Parent'] == 'mother'], 'برای بررسی در مادر'),
            (for_check[for_check['Parent'] == 'child'], 'برای بررسی در فرزند'),
            (not_shared_path, 'موارد پاتوژن غیرمشترک در فرزند و پدر و مادر')
        ]

        self.save_xlsx(datasets, self.output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--keep-intronic', action='store_true')
    parser.add_argument('--no-keep-intronic', dest='keep-intronic', action='store_false')
    parser.set_defaults(keep_intronic=False)

    subparsers = parser.add_subparsers(dest='mode', required=True)

    father_mother = subparsers.add_parser('father_mother', help="Filter the father_mother dataset")
    father_mother.add_argument('mother', help="The file address for the mother gene csv file")
    father_mother.add_argument('father', help="The file address for the father gene csv file")
    father_mother.add_argument('mother_path', help="The file address for the mother pathogen csv file")
    father_mother.add_argument('father_path', help="The file address for the father pathogen csv file")
    
    mother_child = subparsers.add_parser('mother_child', help="Filter the mother_child dataset")
    mother_child.add_argument('mother', help="The file address for the mother gene csv file")
    mother_child.add_argument('child', help="The file address for the child gene csv file")
    mother_child.add_argument('mother_path', help="The file address for the mother pathogen csv file")
    mother_child.add_argument('child_path', help="The file address for the child pathogen csv file")
    mother_child.add_argument('omim', help="The file address for the omim txt file")
    
    father_mother_child = subparsers.add_parser('father_mother_child', help="Filter the father_mother_child dataset")
    father_mother_child.add_argument('father', help="The file address for the father gene csv file")
    father_mother_child.add_argument('mother', help="The file address for the mother gene csv file")
    father_mother_child.add_argument('child', help="The file address for the child gene csv file")
    father_mother_child.add_argument('father_path', help="The file address for the father pathogen csv file")
    father_mother_child.add_argument('mother_path', help="The file address for the mother pathogen csv file")
    father_mother_child.add_argument('child_path', help="The file address for the child pathogen csv file")
    father_mother_child.add_argument('omim', help="The file address for the omim txt file")

    args = parser.parse_args()
    
    if args.mode == 'father_mother':
        file_name = generate_file_name(args.father, args.mother)
        FatherMotherParser(args.mother, args.father, args.mother_path, args.father_path, file_name, args.keep_intronic).run()
    elif args.mode == 'mother_child':
        file_name = generate_file_name(args.child, args.mother)
        MotherChildParser(args.mother, args.child, args.mother_path, args.child_path, args.omim, file_name, args.keep_intronic).run()
    elif args.mode == 'father_mother_child': 
        file_name = generate_file_name(args.father, args.child, args.mother)
        FatherMotherChildParser(args.mother, args.father, args.child, args.mother_path, args.father_path, args.child_path, args.omim, file_name, args.keep_intronic).run()

if __name__ == "__main__":
    main()
