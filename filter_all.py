import pandas as pd
import argparse
import io
import pathlib

from typing import Callable, Tuple, List, Union

Path = Union[str, pathlib.Path]

def generate_file_name(*files: List[str]) -> str: 
    files = [file.split('/')[-1].replace('filtered_', '').replace('.csv','')[:-2] for file in files]
    return f"filtered_{'&'.join(files)}.xlsx"

class GeneralParser():
    match_columns = ["Het Iranome", "Hom Iranome", "Het Our DB", "Chr", "Start", "End", "Ref", "Alt", "Zygosity", "Gene.refGene", "ExonicFunc.refGene"]
    gene_exceptions = ['frameshift insertion', 'frameshift deletion', 'stopgain', 'stoploss', 'splicing']
    
    def warn(self, message: str):
        print("Warning:", message)

    def error(self, message: str):
        print("Error:", message)

    def save_xlsx(self, dataframes: List[Tuple[pd.DataFrame, str]], output: Path):
        current_row = 1
        writer = pd.ExcelWriter(output, engine='xlsxwriter')

        # Write the first dataframe to create the Sheet
        first_df, first_label = dataframes[0]
        first_df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=current_row, startcol=0)
        current_row += first_df.shape[0]

        del dataframes[0]

        # Get the worksheet and create the style
        worksheet = writer.sheets['Sheet1']
        workbook = writer.book

        merge_format = workbook.add_format({
            'bold':     True,
            'align':    'center',
            'valign':   'vcenter',
            'fg_color': 'black',
            'font_color': 'white',
            'font_size': 16
        })

        worksheet.merge_range(0, 0, 0, 3, first_label, merge_format)
        current_row += 1

        # Write the rest of the dataframes
        for (df, label) in dataframes:
            worksheet.merge_range(current_row, 0, current_row, 3, label, merge_format)
            current_row += 1
            df.to_excel(writer, sheet_name='Sheet1', header=False, index=False, startrow=current_row)
            current_row += df.shape[0]
            current_row += 1

        writer.save()

    def read_csv(self, path: Path, parent: str, data_filter: Callable[[pd.DataFrame, str], pd.DataFrame]) -> pd.DataFrame:
        try:
            df = pd.read_csv(path)
        except:
            self.warn('Warning: File has contains bad lines, Some data might get corrupted!')
            df = pd.read_csv(path, engine='python', on_bad_lines=lambda x: x[:273])
            df = df.shift(2, axis='columns').fillna('.')

        return data_filter(df, parent)

    def concat_dataframes(self, dataframes: List[Tuple[pd.DataFrame, str]]):
        dataframes = [df.astype(str).reset_index(drop=True) for df in dataframes]
        return pd.concat(dataframes).sort_values(['Gene.refGene'])

    def mother_and_father_share_gene(self, row: pd.Series):
        row = row.sort_values('Parent')
        parents = tuple(row['Parent'].values)
        return ('mother' in parents) and ('father' in parents)

    def mother_and_child_share_gene(self, row: pd.Series):
        row = row.sort_values('Parent')
        return (tuple(row['Parent'].values) == ('child', 'mother')) and (tuple(row['Zygosity'].values) == ('hom', 'het'))

    def mother_father_and_child_do_not_share_gene(self, row: pd.Series):
        row = row.sort_values('Parent')
        return ('father' not in row['Parent'].values) and ('mother' not in row['Parent'].values)

    def mother_and_child_or_father_and_child_share_gene(self, row : pd.Series):
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

    def filter_normal(self, df: pd.DataFrame, parent: str) -> pd.DataFrame:
        df = df[(df['Hom Iranome'] == '0') | (df['Hom Iranome'] == '.')]
        
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

    def filter_path(self, df: pd.DataFrame, parent: str) -> pd.DataFrame:
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
                       output: Path):
        self.mother = mother
        self.father = father
        self.mother_path = mother_path
        self.father_path = father_path
        self.output = output
    
    def run(self): 
        mother = self.read_csv(self.mother, 'mother', self.filter_normal)
        father = self.read_csv(self.father, 'father', self.filter_normal)

        mother_path = self.read_csv(self.mother_path, 'mother', self.filter_path)
        father_path = self.read_csv(self.father_path, 'father', self.filter_path)

        normal_df = self.concat_dataframes([mother, father])
        path_df = self.concat_dataframes([mother_path, father_path])

        mother_and_father_shared_gene = normal_df.groupby(self.match_columns).filter(self.mother_and_father_share_gene)
        compound_gene = normal_df[(normal_df.duplicated(subset=['Gene.refGene'], keep=False))]
        dangerous_gene = normal_df[(normal_df['ExonicFunc.refGene'].isin(self.gene_exceptions))
                          | (normal_df['ExonicFunc.ensGene'].isin(self.gene_exceptions))
                          | (normal_df['ExonicFunc.knownGene'].isin(self.gene_exceptions))
                          | (normal_df['Func.refGene'].isin(self.gene_exceptions))
                          | (normal_df['Function_description'].isin(self.gene_exceptions))]
        
        mother_and_father_shared_path = path_df.groupby(self.match_columns).filter(self.mother_and_father_share_gene)
        mother_and_father_not_shared_path = path_df[~(path_df.duplicated(subset=self.match_columns, keep=False))]
        
        mother_hom_check = path_df[(path_df['Zygosity'] == 'hom') & (path_df['Parent'] == 'mother')]
        father_hom_check = path_df[(path_df['Zygosity'] == 'hom') & (path_df['Parent'] == 'father')]

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
                       output: Path):
        self.mother = mother
        self.child = child
        self.mother_path = mother_path
        self.child_path = child_path
        self.output = output
    
    def run(self): 
        mother = self.read_csv(self.mother, 'mother', self.filter_normal)
        child = self.read_csv(self.child, 'child', self.filter_normal)

        mother_path = self.read_csv(self.mother_path, 'mother', self.filter_path)
        child_path = self.read_csv(self.child_path, 'child', self.filter_path)

        normal_df = self.concat_dataframes([mother, child])
        path_df = self.concat_dataframes([mother_path, child_path])

        shared_mother_child_gene = normal_df.groupby(self.match_columns).filter(self.mother_and_child_share_gene)
        shared_mother_child_path = path_df.groupby(self.match_columns).filter(self.mother_and_child_share_gene)
        not_shared_mother_child_gene = normal_df.groupby('Gene.refGene').filter(mother_father_and_child_do_not_share_gene)
        not_shared_mother_child_path = path_df.groupby('Gene.refGene').filter(mother_father_and_child_do_not_share_gene)
        
        shared_mother_child = self.concat_dataframes([
            shared_mother_child_gene, 
            shared_mother_child_path
        ])

        not_shared_mother_child = self.concat_dataframes([
            not_shared_mother_child_gene,
            not_shared_mother_child_path
        ])

        dataframes = [
            (shared_mother_child, 'موارد مشترک در مادر و فرزند'),
            (not_shared_mother_child, 'موارد غیرمشترک در فرزند (‌فقط فرزند)')
        ]

        self.save_xlsx(dataframes, self.output)



class FatherMotherChildParser(GeneralParser): 
    def __init__(self, 
                 mother: Path, father: Path, child: Path,
                 mother_path: Path, father_path: Path, child_path: Path,
                 output: Path): 
        self.mother = mother
        self.father = father
        self.child = child
        self.mother_path = mother_path
        self.father_path = father_path
        self.child_path = child_path
        self.output = output

    def run(self): 
        father = self.read_csv(self.father, 'father', self.filter_normal)
        mother = self.read_csv(self.mother, 'mother', self.filter_normal)
        child = self.read_csv(self.child, 'child', self.filter_normal)

        father_path = self.read_csv(self.father_path, 'father', self.filter_path)
        mother_path = self.read_csv(self.mother_path, 'mother', self.filter_path)
        child_path = self.read_csv(self.child_path, 'child', self.filter_path)

        normal_df = self.concat_dataframes([father, mother, child])
        path_df = self.concat_dataframes([father_path, mother_path, child_path])

        shared_gene = normal_df.groupby('Gene.refGene').filter(self.mother_and_child_or_father_and_child_share_gene)
        shared_gene_path = path_df.groupby('Gene.refGene').filter(self.mother_and_child_or_father_and_child_share_gene)
        
        not_shared_gene = normal_df.groupby('Gene.refGene').filter(self.mother_father_and_child_do_not_share_gene)
        not_shared_gene_path = path_df.groupby('Gene.refGene').filter(self.mother_father_and_child_do_not_share_gene)

        datasets = [
            (self.concat_dataframes([shared_gene, shared_gene_path]), 'موارد مشترک در پدر و مادر و فرزند'),
            (self.concat_dataframes([not_shared_gene, not_shared_gene_path]), 'موارد مشترک در پدر و مادر و فرزند')
        ]

        self.save_xlsx(datasets, self.output)


def main():
    parser = argparse.ArgumentParser()
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
    
    father_mother_child = subparsers.add_parser('father_mother_child', help="Filter the father_mother_child dataset")
    father_mother_child.add_argument('father', help="The file address for the father gene csv file")
    father_mother_child.add_argument('mother', help="The file address for the mother gene csv file")
    father_mother_child.add_argument('child', help="The file address for the child gene csv file")
    father_mother_child.add_argument('father_path', help="The file address for the father pathogen csv file")
    father_mother_child.add_argument('mother_path', help="The file address for the mother pathogen csv file")
    father_mother_child.add_argument('child_path', help="The file address for the child pathogen csv file")
    
    args = parser.parse_args()
    
    if args.mode == 'father_mother':
        file_name = generate_file_name(args.father, args.mother)
        FatherMotherParser(args.mother, args.father, args.mother_path, args.father_path, file_name).run()
    elif args.mode == 'mother_child':
        file_name = generate_file_name(args.child, args.mother)
        MotherChildParser(args.mother, args.child, args.mother_path, args.child_path, file_name).run()
    elif args.mode == 'father_mother_child': 
        file_name = generate_file_name(args.father, args.child, args.mother)
        FatherMotherChildParser(args.mother, args.father, args.child, args.mother_path, args.father_path, args.child_path, file_name).run()

if __name__ == "__main__":
    main()