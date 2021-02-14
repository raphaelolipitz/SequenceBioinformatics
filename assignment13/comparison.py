from typing import Dict, List
import pandas as pd

def parse_silva_csv(path: str, header = True, min_len = 6) -> Dict[str, List[str]]:
    with open(path) as f:
        d = {}
        if (header):
            f.readline()
        for line in f:
            ls = []
            s = ''
            incell = False
            for c in line:
                if c == ';': #if this is a separating ; save s and start new cell
                    if (not incell):
                        if (s[-1] == ';'):
                            s = s[:-1]
                        ls.append(s)
                        s = ''
                        continue
                if c == '"': #A " might confine cells
                    incell = not incell
                    continue
                if (c == '\n' and s != ''):
                    if (s[-1] == ';'):
                        s = s[:-1]
                    ls.append(s)
                    break
                s += c
                
            classif = ls[1].split(';')
            while len(classif) < min_len:
                classif.append('Unclassified')
            d[ls[0]] = classif
    return d

def parse_RDP_txt(path: str, min_len = 6) -> Dict[str, List[str]]:
    with open(path) as f:
        d = {}
        in_header = True
        for line in f:
            # skip header, which is separated with a blank line
            if line[0] == '\n':
                in_header = False
                continue
            
            if (not in_header):
                
                ls = []
                s = ''
                for c in line:
                    if c == ';': 
                        ls.append(s)
                        s = ''
                        continue
                    
                    if (c == '\n' and s != ''):
                        ls.append(s)
                        break
                    s += c
                    
                # post-process into the dictionary
                rest = []
                if (len(ls) % 2 == 1): # uneven length -> something went wrong
                    raise RuntimeError
                else:
                    for i in range(int(len(ls) / 2)):
                        if ls[i * 2 + 1] == '+': # + separates ID from the rest
                            key = ls[i * 2]
                            continue
                        
                        if ls[i * 2 + 1][-1] == '%': # this is a confidence
                            if int(ls[i * 2 + 1][:-1]) >= 80:
                                rest.append(ls[i * 2])
                        
                rest = rest[1:] # remove 'Root'
                while len(rest) < min_len:
                    rest.append('Unclassified')
                d[key] = rest
    return d

def compare_names(name1, name2):
    synonyms = [('Bacteroidetes', 'Bacteroidota'),
                ('Actinobacteria', 'Actinobacteriota'),
                ('Verrucomicrobia', 'Verrucomicrobiota')]
    if (name1 == name2):
        return True
    
    for tup in synonyms:
        if name1 in tup and name2 in tup:
            return True
    
    else:
        return False
             
def produce_comparison(RDP, Silva, tax_level):
    ls =[[key, RDP[key][tax_level], 
          Silva[key][tax_level], 
                       compare_names(RDP[key][tax_level], Silva[key][tax_level])]
                       for key in RDP]
        
    df = pd.DataFrame(ls , columns = ['read', 'RDP_class',
                                                   'Silva_class', 'agreement'])
    return df
    
path = 'C:/Users/emilp/Documents/Uni/Sequence Bioinformatics/Assignments/SequenceBioinformatics/assignment13/'

path_32_silva = 'reads16S-32.Silva.csv'
path_500_silva = 'reads16S-500.Silva.csv'
path_32_RDP = 'reads16S-32.RDP.allrank.txt'
path_500_RDP = 'reads16S-500.RDP.allrank.txt'

dRDP32 = parse_RDP_txt(path + path_32_RDP)
dRDP500 = parse_RDP_txt(path + path_500_RDP)
dSilva32 = parse_silva_csv(path + path_32_silva)
dSilva500 = parse_silva_csv(path + path_500_silva)

d32 = (dRDP32, dSilva32)
d500 = (dRDP500, dSilva500)

phylum = 1
order = 3
genus = 5

#df = produce_comparison(dRDP500, dSilva500, order)

result = pd.DataFrame(index = [32,500], columns=['phylum', 'order', 'genus'])

for dataset in (d32, d500):
    for rank in ((phylum, 'phylum'), (order, 'order'), (genus, 'genus')):
        df = produce_comparison(dataset[0], dataset[1], rank[0])
        result.loc[len(df), rank[1]] = len(df.loc[df['agreement'],])
        print('The systems agree on {} of {} reads regarding {}.'
              .format(len(df.loc[df['agreement'],]), len(df), rank[1]))
        
print(result.to_latex())
print()
print(produce_comparison(dRDP32, dSilva32, phylum).to_latex(index = False))











