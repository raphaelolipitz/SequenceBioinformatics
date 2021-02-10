from typing import Dict, List

def parse_silva_csv(path: str) -> Dict[str, List[str]]:
    with open(path) as f:
        d = {}
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
                
            d[ls[0]] = ls[1].split(';')
    return d

def parse_RDP_txt(path: str) -> Dict[str, List[str]]:
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
                        
                d[key] = rest
    return d
                    
path = 'C:/Users/emilp/Documents/Uni/Sequence Bioinformatics/Assignments/SequenceBioinformatics/assignment13/'
path_32_silva = 'reads16S-32.Silva.csv'
path_500_silva = 'reads16S-500.Silva.csv'
path_32_RDP = 'reads16S-32.RDP.allrank.txt'
path_500_RDP = 'reads16S-500.RDP.allrank.txt'

dRDP32 = parse_RDP_txt(path + path_32_RDP)
dRDP500 = parse_RDP_txt(path + path_500_RDP)
dSilva32 = parse_silva_csv(path + path_32_silva)
dSilva500 = parse_silva_csv(path + path_500_silva)

print(dSilva32)
