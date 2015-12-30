//A couple functions for converting numbers to strings

string IntToStr(int n)
{
    stringstream result;
    result << n;
    return result.str();
}

string DoubleToStr(double n)
{
    stringstream result;
    result << n;
    return result.str();
}