struct Vertex
{
        double x;
};
struct Element
{
        struct Vertex *n;
	struct Vertex *pot;
        double h;
};
