__author__ = 'chenzhi'



def extract(src_path, dest_path, digit):
    """
    extract data of specific digit from source data
    :param src_path: path of the source file
    :param dest_path: path of destination file
    :param digit: the digit you want to extract
    :return:
    """
    content = ""
    dest_file = open(dest_path, "w")
    with open(src_path) as f:
        content = f.readlines()
    for i in range(len(content)):
        if content[i] == " " + str(digit) + "\n":
            for j in range(i-32, i):
                dest_file.write(content[j])
            dest_file.write("\n")
    dest_file.close()


if __name__ == "__main__":
    extract("/home/chenzhi/data/optdigits-orig.tra","/home/chenzhi/data/3", 3)