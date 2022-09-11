import sys, os, glob


def main():
    wdir = sys.argv[1]
    head = True
    for fpath in glob.glob(os.path.join(wdir, "*", "ASGAL", "events.wpsi.csv")):
        fcontent = open(fpath).readlines()
        if head:
            print(fcontent[0], end="")
            head = False
        print(fcontent[1], end="")


if __name__ == "__main__":
    main()
