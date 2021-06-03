using System;
using System.Linq;
using Xunit;

namespace DelaunatorNet.Test
{
    public class DelaunatorTest
    {
        [Fact]
        public void DelaunatorBuildTest_ShouldSucceed()
        {
            var points = new double[]{4, 1,
                3.7974166882130675, 2.0837249985614585,
                3.2170267516619773, 3.0210869309396715,
                2.337215067329615, 3.685489874065187,
                1.276805078389906, 3.9872025288851036,
                0.17901102978375127, 3.885476929518457,
                -0.8079039091377689, 3.3940516818407187,
                -1.550651407188842, 2.5792964886320684,
                -1.9489192990517052, 1.5512485534497125,
                -1.9489192990517057, 0.44875144655029087,
                -1.5506514071888438, -0.5792964886320653,
                -0.8079039091377715, -1.394051681840717,
                0.17901102978374794, -1.8854769295184561,
                1.276805078389902, -1.987202528885104,
                2.337215067329611, -1.6854898740651891,
                3.217026751661974, -1.021086930939675,
                3.7974166882130653, -0.08372499856146409};

            Delaunator delaunator = new Delaunator(points);
            var result = delaunator.Build();

            var expectedTriangles = new[] {
                16, 2, 0, 0, 2, 1, 10, 2,16,
                16, 12, 10, 10, 3, 2, 15, 12,
                16, 5, 4, 3, 14, 12, 15, 10,
                5, 3, 13, 12, 14, 7, 6, 5, 12,
                11, 10, 10, 7, 5, 10, 8, 7, 10, 9, 8};
            var expectedHalfEdeges = new int[] { 7, 3, -1, 1, -1, -1, 14, 0, 11, 16, 35, 
                8, 26, -1, 6, 22, 9, -1, -1, -1, 25, 28, 15, -1,
                38, 20, 12, -1, 21, -1, -1, -1, 37, 
                -1, -1, 10, 41, 32, 24, 44, -1, 36, -1, -1, 39 };

            var expectedHull = new int[] { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 16, 15, 14, 13, 12, 11 };

            Assert.Equal(expectedTriangles, result.Triangles);
            Assert.Equal(expectedHull, result.ConvexHull);
            Assert.Equal(expectedHalfEdeges, result.HalfEdges);
            Assert.Equal(points, result.Points);
        }
    }
}
