from defs import *
import gzip
import tarfile
import shutil


def targz_dir(outer_dir, dirs_list, dest_zip_filename, delete_after_zipping):
    """
    Compress data from multiple directories to a single .targz file
    :param outer_dir: where the output file should be written to
    :param dirs_list: the directories names to be compressed
    :param dest_zip_filename: name of output file
    :param delete_after_zipping: True or False
    :return: NA
    """
    cwd = os.getcwd()
    os.chdir(outer_dir)
    tarw = tarfile.open(dest_zip_filename, "w:gz")
    for dirname in dirs_list:
        if os.path.exists(dirname):
            tarw.add(dirname)
    tarw.close()
    if delete_after_zipping:
        for dirname in dirs_list:
            try:
                shutil.rmtree(dirname)
            except Exception as e1:
                print(e1)
                pass
    os.chdir(cwd)


def untargz(zip_file_dest, delete_after_extracting=False):
    """
    De-compress file
    :param zip_file_dest: where to unpack the zipped file
    :param delete_after_extracting: True or False
    :return: NA
    """
    dirpath, zip_filename = os.path.split(zip_file_dest)
    cwd = os.getcwd()
    os.chdir(dirpath)
    tarx = tarfile.open(zip_file_dest, "r:gz")
    tarx.extractall(dirpath)
    tarx.close()
    os.chdir(cwd)
    if delete_after_extracting:
        os.remove(zip_file_dest)


def average(lst):
    """
    Calculates the average of a given list
    :param lst: list of numbers
    :return: average value
    """
    return sum(lst) / len(lst)


def range_of_lst(lst):
    """
    calculates the range of a given list
    :param lst: list of numbers
    :return: max-min
    """
    return max(lst) - min(lst)


def fix_simulated_tree_file(tree_file):
    """
    a patch that adds a semicolon at the end of .phr tree (simulated tree) so the can tree can be read
    :param tree_file:
    :return: NA
    """
    with open(tree_file, "a") as add:
        add.write(";")


def create_counts_hash(counts_file):
    """
    puts all counts from a counts file in a dictionary: taxa_name: count. Taxa with an X count are skipped.
    If there are counts that are X the function returns two hashes and a set of the taxa to prune.
    :param counts_file in FASTA format
    :return:(2) dictionary of counts
    """
    d = {}
    with open(counts_file, "r") as counts_handler:
        for line in counts_handler:
            line = line.strip()
            if line.startswith('>'):  # taxon name
                name = line[1:]
            else:
                if line != "x":
                    num = int(line)
                    d[name] = num
    return d


def regex_internal(str):
    """
    extracts the number from the internal node's name from a phylogeny, in a NX-XX format. Used in get_max_transition(tree_file)
    :param str: internal node's label
    :return: count of node
    """
    tmp = re.search("N\d+\-(\d+)", str)
    if tmp:
        num = int(tmp.group(1))
    return num


def regex_tip(str):
    """
    extracts the number from a tip label in a phylogeny
    :param str: tip label
    :return: count of tip
    """
    tmp = re.search("(\d+)", str)
    if tmp:  # there is a number at the tip, and not X
        num = int(tmp.group(1))
    return num


def get_max_transition(tree_file):
    """
    searches for the largest transition that was made on the phylogeny itself, between internal node and tip, or between internal nodes. Used in base_num_models.py
    :param tree_file: phylogeny file
    :return: a number representing the maximal transition on the tree
    """
    t = Tree(tree_file, format=1)
    max_transition = 0
    for node in t.traverse():
        if node.name == "":
            continue
        if not node.is_leaf():
            num1 = regex_internal(node.name)
            for child in node.get_children():
                if child.is_leaf():  # if the child is a tip - parse number from tip label
                    num2 = regex_tip(child.name)
                else:  # if the child is an internal node - take number
                    num2 = regex_internal(child.name)
                tmp_score = abs(num1 - num2)
                if max_transition < tmp_score:
                    max_transition = tmp_score
    return max_transition


def copy_files(src, dest, files=None):
    """
    copy all files from source directory to destination directory. Used when re-running BASE NUM models, to keep previous results
    :param src: source directory
    :param dest: destination directory
    :param files: a list of specific files to be copied
    :return: NA
    """
    src_files = os.listdir(src)
    if not os.path.exists(dest):
        os.system("mkdir " + dest)
    if files is not None:
        src_files = files
    for file_name in src_files:
        full_src_file_name = src + file_name
        full_dest_file_name = dest + file_name
        if os.path.isfile(full_src_file_name):
            shutil.copy(full_src_file_name, full_dest_file_name)


def extract_line_from_file(filename, str_search, num=False, integer=False):
    """
    uses regular expression to search a string in a file
    :param filename: the file to be searched in
    :param str_search: the string to search for
    :param num: should a number be extracted?
    :param integer: should an integer be returned? if False and num=True then a float will be returned
    :return: True/False for a sentence
                integer if num=True and integer=True
                float if num=True and integer=False
    """
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if num:
                tmp = re.search(str_search + ".*?(\d+)", line)
                if tmp:
                    return int(tmp.group(1)) if integer else float(tmp.group(1))
            else:
                tmp = re.search(str_search, line)
                if tmp:
                    return True
        return False


def n_folders_lst(n):
    """
    used to targz all simulations to a single file
    :return: list with names of folders from 0 to n
    """
    lst = []
    for i in range(n):
        lst.append(str(i))
    return lst


def print_error(msg, exception="", sep="*****"):
    print(sep)
    print(msg, "\n")
    print(sep)
    if exception is not "":
        print(exception)
        print(sep)
