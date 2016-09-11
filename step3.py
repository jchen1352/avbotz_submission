n=0;i=open('test.in.txt');o=open('test.out.txt','w');j=len(i.readlines());i.seek(0);exec"a=int(sum([int(x)**2 for x in i.readline().split()])**.5);n+=a>200;o.write(`a`+'\\n');"*j;o.write(`n`)
