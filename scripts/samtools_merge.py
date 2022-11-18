import sys
import getopt
import subprocess
import os



def chunker_list(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def main(argv):
    print("script")
    inputfile = ''
    outputfile = ''
    n_files = 1000
    try:
        opts, args = getopt.getopt(argv, "hi:o:",
                                   ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('samtools_merge.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('samtools_merge.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    outputfiletmp = f"{outputfile}.tmp"
    chunks = []
    with open(inputfile, "r") as f:
        lines = f.read().splitlines() 

        chunks = list(chunker_list(lines, n_files))
    #print(chunks)
    #print(len(chunks))

    count = 0
    if os.path.exists(outputfile):
        os.remove(outputfile)
    # if not os.path.exists(f"{outputfile}tmp"):
    #     os.makedirs(f"{outputfile}_tmp")

    for chunk in chunks:
        if count == 0:
            files = " ".join(chunk)
            com = f"samtools merge -o {outputfile} {files}"
            #print(com)
            exe = subprocess.call(com, shell=True)
            os.rename(outputfile, outputfiletmp)
        else:
            files = " ".join(chunk)
            com = f"samtools merge -o {outputfile} -f {outputfiletmp} {files}"
            #print(com)
            exe = subprocess.call(com, shell=True)
            os.rename(outputfile, outputfiletmp)
        count += 1
    os.rename(outputfiletmp, outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])