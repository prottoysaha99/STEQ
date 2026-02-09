import sys


def quartetScore(filename):
    id = 0
    with open(filename,"r") as fi:
        for ln in fi:
            if ln.startswith("Final quartet score is: "):
                id=int(ln.split()[4])
                break
    
    return id


def main():
    filepath = str(sys.argv[1])
    itr =int(sys.argv[2])
    max_score = 0
    fileid = 0

    for i in range(0,itr+1):
        filename=filepath + "logs." + str(i)
        score = quartetScore(filename)

        if score > max_score:
            max_score = score
            fileid = i

    return fileid

        

fid = main()
print(fid)
