import sys

lines = [str(n.strip()) for n in sys.stdin]

head = lines[0].split('\t')
head_chr = head[0]
head_start = int(head[1])
head_stop = int(head[2])

for lx in lines[1:]:
    lx = lx.split('\t')
    next_chr = lx[0]
    next_start = int(lx[1])
    next_stop = int(lx[2])
    if head_chr == next_chr:
        if (min([next_start,next_stop]) <= head_start <= max([next_start,next_stop])) or (min([next_start,next_stop]) <= head_stop <= max([next_start,next_stop])):
            head_start = min([head_start,head_stop,next_start,next_stop])
            head_stop = max([head_start,head_stop,next_start,next_stop])
        else:
            print('\t'.join(map(str,[head_chr,head_start,head_stop])))
            head_chr = next_chr 
            head_start = next_start
            head_stop = next_stop
    else:
        print('\t'.join(map(str,[head_chr,head_start,head_stop])))
        head_chr = next_chr
        head_start = next_start
        head_stop = next_stop

print('\t'.join(map(str,[head_chr,head_start,head_stop])))

