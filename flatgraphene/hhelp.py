from pathlib import Path

def help():
    """
    Forward documentation/doc strings to the terminal
    """
    #construct path to help document from current file
    help_file_name = 'help_doc.md'
    help_file_path = str(Path(__file__).parent) + '/' + help_file_name

    #read lines of help file
    with open(help_file_path,'r') as f:
        lines = f.readlines()

    #send lines of help file to terminal
    for line in lines:
        print(line.strip()) #remove new line character as print adds its own

    
