import sys, string, os

if not os.path.exists("logs"):
    os.mkdir("logs")

for i in range(3,15):
    with open("data", "rt") as fin:
        with open("data_", "wt") as fout:
            for line in fin:
                namefile = "tree" + str(i) + ".xml"
                fout.write(line.replace('\txmlfile = tree3.xml', '\txmlfile = ' + namefile))

        os.system("./NavierStokesExample data_ > logs/log" + str(i) + ".txt")
        os.remove("data_")
