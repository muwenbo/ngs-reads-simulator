import os


class ReadsFileHandle:
    """
    A class to handle file operations for simulated reads and associated data.

    This class manages the creation, opening, closing, and deletion of files
    used for storing simulated reads, variant information, and error data.
    """

    def __init__(self, file_prefix, var_out_file, error_out_file):
        """
        Initialize the ReadsFileHandle object.

        Args:
            file_prefix (str): Prefix for the read files.
            var_out_file (str): Name of the file to store variant information.
            error_out_file (str): Name of the file to store error information.
        """
        self.file_names = [f"{file_prefix}1.fastq", f"{file_prefix}2.fastq"]
        self.var_out_file = var_out_file
        self.error_out_file = error_out_file

    def file_open(self):
        """
        Open all necessary files for writing.
        """
        self.error_fh = open(self.error_out_file, 'w')
        self.var_fh = open(self.var_out_file, 'w')
        self.reads1_f = open(self.file_names[0], 'a')
        self.reads2_f = open(self.file_names[1], 'a')
            
    def file_close(self):
        """
        Close all opened files.
        """
        self.error_fh.close()
        self.var_fh.close()
        self.reads1_f.close()
        self.reads2_f.close()
    
    def file_delete(self):
        """
        Delete all associated files if they exist.
        """
        for file in self.file_names:
            if os.path.isfile(file):
                print(f"Removing {file}!")
                os.remove(file)
        
        if os.path.isfile(self.var_out_file):
            print(f"Removing {self.var_out_file}!")
            os.remove(self.var_out_file)

        if os.path.isfile(self.error_out_file):
            print(f"Removing {self.error_out_file}!")
            os.remove(self.error_out_file)

