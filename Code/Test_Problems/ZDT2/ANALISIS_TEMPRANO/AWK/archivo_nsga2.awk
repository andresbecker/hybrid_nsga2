BEGIN{getline eval < "eval.out"}
{
    if(NR==1)
    {
        popsize=$0;
        print $0;
    }
    if(NR==2)
        printf("%d\n",eval/popsize);
    if(NR==3)
        print $0;
    if(NR==4)
        print $0;
    if(NR==5)
        final=5+$0+14;
    if(NR>2 && NR<final)
        print $0;
    if(NR==final)
        printf("%d\n",eval/popsize + 2);
}
