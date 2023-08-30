
static double f1f2(double x1, double x2)
{
    return Math.Round((-2 * Math.Pow(x1, 3) -4 * Math.Pow(x1, 2) -24 * x1 + 17) * (2 * Math.Pow(x2, 2) - 10 * x2 + 15), 8);
}

static double f2f1(double x1, double x2)
{
    return Math.Round((3 * Math.Pow(x2, 2) -18 * x2 * x1 -33 * x2 - 12) * (Math.Pow(x1, 2) - 6 * x1 + 13), 8);
}


double x1LowerEnd = 1;
double x1HigherEnd = 8;
double x2LowerEnd = 0;
double x2HigherEnd = 6;

const double x1Step = 0.02;
const double x2Step = 0.02;

List<Table> tables = new();
Console.WriteLine(String.Format("|{0, 8}|{1, 8}|{2, 20}|{3, 20}|", "x1", "x2", "f12", "f21"));
for (double x1 = x1LowerEnd; x1 <= x1HigherEnd; x1 = Math.Round(x1 + x1Step, 8))
{
    for (double x2 = x2LowerEnd; x2 <= x2HigherEnd; x2 = Math.Round(x2 + x2Step, 8))
    {
        Table table = new(x1, x2, f1f2(x1, x2), f2f1(x1, x2));
        tables.Add(table);
        Console.WriteLine(table);
    }
}

Table? guaranteedF12 = tables.GroupBy(q => q.X1).Select(q => q.MinBy(g => g.F12)).MaxBy(q => q.F12);
Table? guaranteedF21 = tables.GroupBy(q => q.X2).Select(q => q.MinBy(g => g.F21)).MaxBy(q => q.F21);

if(guaranteedF12 is not null && guaranteedF21 is not null)
{
    Console.WriteLine($"Guaranteed result f12* = {guaranteedF12.F12}:");
    Console.WriteLine(String.Format("|{0, 8}|{1, 8}|{2, 20}|{3, 20}|", "x1", "x2", "f12", "f21"));
    Console.WriteLine(guaranteedF12);
    Console.WriteLine($"Guaranteed result f21* = {guaranteedF21.F21}:");
    Console.WriteLine(String.Format("|{0, 8}|{1, 8}|{2, 20}|{3, 20}|", "x1", "x2", "f12", "f21"));
    Console.WriteLine(guaranteedF21);
}


class Table
{
    public double X1 { get; set; }
    public double X2 { get; set; }
    public double F12 { get; set; }
    public double F21 { get; set; }

    public Table(double x1, double x2, double f12, double f21)
    {
        X1 = x1;
        X2 = x2;
        F12 = f12;
        F21 = f21;
    }

    public override string ToString()
    {
        return String.Format("|{0, 8}|{1, 8}|{2, 20}|{3, 20}|", X1, X2, F12, F21);
    }
}