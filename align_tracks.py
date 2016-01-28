def main():
    pass

if __name__ == '__main__':
    chr1_file = open('./E116_15_coreMarks_chr1_statebyline.txt')
    # Read in all the lines at once
    contents = chr1_file.read()
    # Dump the first 2 lines since we don't need them
    track = contents.split('\n')[2:]
    print sum(1 for x in track if x == '7')
    bedfile = open('./wgEncodeBroadHistoneGm12878H3k27acStdSig.bedGraph')
    
    # enhancers = []
    # curr = []
    # next_point = 0
    # while True:
        # Find the next group of 7s
        # while track[next_point] != '7':
            # next_point += 1
        # start = next_point
        # while track[next_point] == '7':
            # next_point += 1
        

    enhancers = []
    curr = [(0,0)]
    for line in bedfile:
        tokens = line.strip().split('\t')
        if track[int(tokens[1]) / 200] == '7' or track[(int(tokens[1]) - 100) / 200] == '7' or  track[(int(tokens[1]) + 100) / 200] == '7':
            curr.append((float(tokens[3]), int(tokens[2]) - int(tokens[1]) + curr[-1][1]))
        else:
            if len(curr) > 1:
                enhancers.append(curr)
                curr = [(0,0)]

    print enhancers[:10]
