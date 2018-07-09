import sys
import math

# ______________________________________________________________________
def main ():
    try:
        # assign command line arguments to variables
        filename    = sys.argv[1]
        cutoff      = float(sys.argv[2])
        dcutoff     = float(sys.argv[3])
    except IndexError:
        print ("Incorrect number of arguments (expecting 2)")
        print ("Syntax: python make_dotfile [edgefile] [cutoff] [directional cutoff]")
        sys.exit()
    
    # read contents of file
    infile = open(filename, 'U')
    lines = infile.readlines()
    infile.close()
    
    # prepare DOT file
    outfile = open(filename+'.dot', 'w')
    outfile.write('digraph foo\n{\n\tnode [shape="box" fontname="Helvetica"];\n\trankdir=LR;\n')
    outfile.write('\tedge [labelfontname="Helvetica" labelangle=30 labeldistance=2];\n')    
    
    # set up matrix
    nvar = math.sqrt(len(lines))
    if nvar % 1 > 0:
        print ("ERROR: Number of lines not compatible with square matrix")
        sys.exit()
    
    nvar = int(nvar)
    print ("Detected %d variables in edgelist" % nvar)
    
    # convert edgelist to matrix
    mat = [[0 for row in range(nvar)] for col in range(nvar)]
    labels = []
    for ln in range(len(lines)):
        row = int(ln/nvar)
        col = int(ln%nvar)
        items = lines[ln].strip('\n').split(',')
        
        if items[0] not in labels:
            labels.append(items[0])
        
        mat[row][col] = float(items[2])
    
    #print (labels)
    
    # generate list of nodes involved in significant edges, not including response or race
    nodes = []
    for i in range(nvar-1):
        for j in range(i+1, nvar):
            prob = mat[i][j] + mat[j][i]
            if prob > cutoff:
                if i not in nodes: nodes.append(i)
                if j not in nodes: nodes.append(j)
    
    nodes.sort()
    
    for i in nodes:
        label = labels[i]
        outfile.write('\t"'+label+'" [style="filled" fillcolor="white" fontcolor="black"];\n')
    
    # output to DOT file
    for i in range(nvar-1):
        for j in range(i+1, nvar):
            
            prob = mat[i][j] + mat[j][i]
            
            if prob > cutoff:
                if mat[j][i]/prob >= dcutoff:
                    outfile.write('\t"'+labels[j]+'"->"'+labels[i]+'" [dir="forward" headlabel=%d];\n' % int(prob*100))
                elif mat[i][j]/prob >= dcutoff:
                    outfile.write('\t"'+labels[i]+'"->"'+labels[j]+'" [dir="forward" headlabel=%d];\n' % int(prob*100))
                else:
                    outfile.write('\t"'+labels[i]+'"->"'+labels[j]+'" [dir="none" headlabel=%d];\n' % int(prob*100))
    
    outfile.write('}\n\n');
    outfile.close()
    
    
# ______________________________________________________________________
if __name__ == "__main__":
    main()
