##########################################
# Takes in manual logging accepts file to track all accepted iterations
# Filters a given log file to only keep accepted iterations
#
# (Likely only relevant if the filtered log file logged every iteration)
##########################################

def main():
    with open(f"manualLoggingAccepts.txt", "r") as inp:
        accepts = set(['0'])
        text = inp.read()
        lines = text.split('\n')
        for line in lines[1:-1]:
            accepts.add(line.split(' ')[0])

    with open (f'testSeedbankBasic.every.log', 'r') as tgt:
        with open(f"processAccepts.output", "w") as out:
            text = tgt.read()
            lines = text.split('\n')
            out.write(lines[0] + "\n")

            for line in lines[1:-1]:
                if (line.split('\t')[0] in accepts):
                    out.write(line + "\n")
                # spl = line.split('\t')
                # if (spl[0] in accepts):
                #     out.write('\t'.join(spl) + '\n')

if __name__ == "__main__":
    main()