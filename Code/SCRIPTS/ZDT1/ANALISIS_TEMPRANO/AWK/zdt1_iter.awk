BEGIN{getline iter < "temp.out"}
{
    if(NR==1)
    {
        print $0;
    }
    if(NR==2)
    {
        print $0;
    }
    if(NR==3)
        print $0;
    if(NR==4)
        print $0;
    if(NR==5)
        final=5+$0+14;
    if(NR>2 && NR<final)
        print $0;
    if(NR==final)
        print $0+iter-1;
}
