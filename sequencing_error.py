class SequencingError:
    """
    A class to represent and manage sequencing errors in simulated reads.

    This class stores information about sequencing errors, including their
    positions and the changes in alleles. It provides methods to add new
    errors and write the error information to a file.
    """

    def __init__(self):
        """
        Initialize the Error object with an empty dictionary to store error information.
        """
        self.error_info = {}

    def add_item(self, ref_chr, ref_start, parent_index, local_index, raw_allele, new_allele):
        """
        Add a new sequencing error to the error information dictionary.

        Args:
            ref_chr (str): Reference chromosome.
            ref_start (int): Start position in the reference.
            parent_index (int): Index of the parent read.
            local_index (int): Local index within the read where the error occurred.
            raw_allele (str): Original allele before the error.
            new_allele (str): New allele after the error.
        """
        global_index = str(ref_start + parent_index + local_index)
        key = f"{ref_chr}:{global_index}"
        
        if key not in self.error_info:
            self.error_info[key] = [ref_chr, global_index, raw_allele, new_allele]
        else:
            # If multiple errors occur at the same position, concatenate the new alleles
            self.error_info[key][3] = f"{self.error_info[key][3]}\\{new_allele}"

    def write_to_file(self, file_h):
        """
        Write all stored error information to a file.

        Args:
            file_h (file): An open file handle to write the error data.
        """
        for key in self.error_info:
            file_h.write('\t'.join(self.error_info[key]) + '\n')
