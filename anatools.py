import pandas as pd
import pickle
import os
import re
import csv
from io import StringIO
import os 

def load_pickle(path,filename='ana.pickle'):
    fullpath = path + '/' + filename
    with open(fullpath,'rb') as file:
        data = pickle.load(file)
    return data

def load_pickles(dfconf,filename='ana.pickle'):
    paths = dfconf['PATH'].values
    datas = []
    for p in paths:
        data = load_pickle(p,filename)
        datas.append(data)
    return datas

def load_analysis(key,dfpickels):
    res = []
    for i,df in enumerate(dfpickels):
        res_ = df[key]
        res.append(res_)
    return res

def load_ana(key,dfpickels,colkeys):
    res = load_analysis(key,dfpickels)
    rez = []
    for r in res:
        r_ = r.loc[:,colkeys]
        rez.append(r_)
    df = pd.concat(rez,axis=1)
    return df
    

def load_rmsd(dfpickels,key='rmsd',colkeys=2):
    return load_ana(key,dfpickels,colkeys)

def load_rg(dfpickels,key='rg',colkeys=2):
    return load_ana(key,dfpickels,colkeys)

def load_sasa(dfpickels,key='sasa',colkeys=2):
    return load_ana(key,dfpickels,colkeys)

def load_hbond(dfpickels,key='hbond',colkeys=2):
    return load_ana(key,dfpickels,colkeys)

def mmpbsa_to_csv(input_filepath: str) -> str:
    """
    Parses the text output of a gmx_MMPBSA calculation from a file and 
    converts it into a CSV file saved in the same directory as the input.

    Args:
        input_filepath: The path to the gmx_MMPBSA output log file (e.g., .dat or .log).

    Returns:
        A confirmation string indicating the output file path, or an error message.
    """
    
    # Read the content from the specified file path
    try:
        with open(input_filepath, 'r') as f:
            mmpbsa_output_text = f.read()
    except FileNotFoundError:
        return f"Error: Input file not found at path: {input_filepath}"
    
    # 1. Define the structure and pattern for extracting sections
    systems = ['Complex', 'Receptor', 'Ligand', 'Delta']
    
    # Use a regex pattern to find the relevant tables for each system.
    data_pattern = re.compile(
        r'(\n|^)(' + '|'.join(systems) + r'|\S+ Delta)[\s\S]*?(?:[-]{60,}\s*?)([\s\S]*?)(?=\n[-]{60,}|\Z)',
        re.MULTILINE | re.IGNORECASE
    )
    
    parsed_data = []
    
    # 2. Extract data for each system
    for match in data_pattern.finditer(mmpbsa_output_text):
        section_prefix = match.group(2).strip()
        table_content = match.group(3).strip()
        
        # Normalize system name
        if 'Delta' in section_prefix:
            system_name = 'Delta'
        else:
            system_name = section_prefix.split(':')[0].strip()

        # Iterate over each line in the table content
        for line in table_content.split('\n'):
            line = line.strip()
            if not line or line.startswith('-'):
                continue

            # Split the line by spaces, filtering out empty strings
            parts = [p for p in line.split() if p]
            
            # The structure is: Component + 5 values
            if len(parts) >= 6:
                component = parts[0]
                try:
                    # Convert values to float for consistency, ensuring correct order
                    values = [float(p) for p in parts[1:6]]
                    
                    parsed_data.append({
                        'System': system_name,
                        'Energy_Component': component,
                        'Average': values[0],
                        'SD(Prop.)': values[1],
                        'SD': values[2],
                        'SEM(Prop.)': values[3],
                        'SEM': values[4],
                    })
                except ValueError:
                    # Skip lines that might look like data but contain non-numeric values
                    continue

    # 3. Write data directly to a CSV file
    
    # Determine output file path by changing the extension to .csv
    base, _ = os.path.splitext(input_filepath)
    output_filepath = base + ".csv"
    
    try:
        # 'newline=' is crucial for proper CSV handling on all operating systems
        with open(output_filepath, 'w', newline='') as outfile:
            fieldnames = ['System', 'Energy_Component', 'Average', 'SD(Prop.)', 'SD', 'SEM(Prop.)', 'SEM']
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    
            writer.writeheader()
            writer.writerows(parsed_data)
        
        return f"Successfully parsed data and saved to CSV file: {output_filepath}"
    
    except Exception as e:
        return f"Error writing CSV file to {output_filepath}: {e}"

def convert_bfedat(bfefname,path,fname='FINAL_RESULTS_MMPBSA'):
    fpath = path + '/' + bfefname +'/' + fname 
    mmpbsa_to_csv(fpath + '.dat')

def load_bfedat(bfefname,path,fname='FINAL_RESULTS_MMPBSA'):
    fpath = path + '/' + bfefname +'/' + fname 
    df = pd.read_csv(fpath + '.csv')
    return df