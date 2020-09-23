BEGIN{getline iter < "temp.out"}
{
    if(NR==2)
    {
        print $0+iter;
    }
    else
    {
	print $0;
    }
}
