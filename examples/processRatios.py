##########################################
# Takes a log file and calculates the difference between each value and itself on the last line
# For log probability values, returns the log ratio. 
#                       --> positive ~ > 1
#                           negative ~ < 1 
# 
##########################################

def main():
    with open(f"processAccepts.output", "r") as inp:
        with open(f"processRatios.Accepts.output", "w") as out:
    # with open(f"testSeedbankBasic.log", "r") as inp:
        # with open(f"processRatios.testSeedbankBasic.output", "w") as out:
            text = inp.read()
            lines = text.split('\n')
            out.write('\t\t\t\t'.join(lines[0].split('\t')) + "\n")
            last_values = lines[1].split('\t')
            for i in range(2, len(lines)-1):
                values = lines[i].split('\t')
                out_values = []
                for j in range(len(values)):
                    out_values.append(str(format(float(values[j]) - float(last_values[j]), '.4f')))
                out_values[0] = values[0]
                out.write("\t\t\t\t\t\t".join(out_values) + "\n")
                last_values = values

if __name__ == "__main__":
    main()