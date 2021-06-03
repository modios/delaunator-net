namespace DelaunatorNet
{
    public class TriangulationInfo
    {
        public int[] ConvexHull { get; }
        public int[] Triangles { get; }
        public double[] Points { get; }
        public int[] HalfEdges { get; }

        public TriangulationInfo(int[] hull, int[] triangles, double[] points, int[] halfEdges)
        {
            ConvexHull = hull;
            Triangles = triangles;
            Points = points;
            HalfEdges = halfEdges;
        }
    }
}
