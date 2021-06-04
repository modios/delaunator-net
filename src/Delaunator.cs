
using RobustPredicates;
using System;

namespace DelaunatorNet
{
    public class Delaunator
    {
        private readonly static double _epsilon = Math.Pow(2, -52);
        private readonly static int[] _edgeStack = new int[512];
        private readonly double[] _coords;
        private readonly int _n;
        private readonly int _maxTriangles;
        private readonly int[] _hullHash;
        private readonly int _hashSize;
        private readonly int[] _hullPrev;
        private readonly int[] _hullNext;
        private readonly int[] _hullTri;
        private readonly int[] _ids;
        private readonly double[] _dists;
        private int _hullStart;
        private int _trianglesLen;
        private int[] _triangles;
        private int[] _halfedges;
        private double _cx;
        private double _cy;

        private int HashKey(double x, double y)
        {
            return (int)Math.Floor(DelaunatorHelpers.PseudoAngle(x - _cx, y - _cy) * _hashSize) % _hashSize;
        }

        private void Link(int a, int b)
        {
            _halfedges[a] = b;
            if (b != -1)
            {
                _halfedges[b] = a;
            }
        }

        private int AddTriangle(int i0, int i1, int i2, int a, int b, int c)
        {
            var t = _trianglesLen;

            _triangles[t] = i0;
            _triangles[t + 1] = i1;
            _triangles[t + 2] = i2;

            Link(t, a);
            Link(t + 1, b);
            Link(t + 2, c);

            _trianglesLen += 3;

            return t;
        }

        private int Legalize(int a)
        {
            int i = 0;
            int ar;

            // recursion eliminated with a fixed-size stack
            while (true)
            {
                var b = _halfedges[a];

                /* if the pair of triangles doesn't satisfy the Delaunay condition
                 * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
                 * then do the same check/flip recursively for the new pair of triangles
                 *
                 *           pl                    pl
                 *          /||\                  /  \
                 *       al/ || \bl            al/    \a
                 *        /  ||  \              /      \
                 *       /  a||b  \    flip    /___ar___\
                 *     p0\   ||   /p1   =>   p0\---bl---/p1
                 *        \  ||  /              \      /
                 *       ar\ || /br             b\    /br
                 *          \||/                  \  /
                 *           pr                    pr
                 */
                var a0 = a - a % 3;
                ar = a0 + (a + 2) % 3;

                if (b == -1)
                { // convex hull edge
                    if (i == 0) break;
                    a = _edgeStack[--i];
                    continue;
                }

                var b0 = b - b % 3;
                var al = a0 + (a + 1) % 3;
                var bl = b0 + (b + 2) % 3;

                var p0 = _triangles[ar];
                var pr = _triangles[a];
                var pl = _triangles[al];
                var p1 = _triangles[bl];

                var illegal = DelaunatorHelpers.InCircle(
                    _coords[2 * p0], _coords[2 * p0 + 1],
                    _coords[2 * pr], _coords[2 * pr + 1],
                    _coords[2 * pl], _coords[2 * pl + 1],
                    _coords[2 * p1], _coords[2 * p1 + 1]);

                if (illegal)
                {
                    _triangles[a] = p1;
                    _triangles[b] = p0;

                    var hbl = _halfedges[bl];

                    // edge swapped on the other side of the hull (rare); fix the halfedge reference
                    if (hbl == -1)
                    {
                        var e = _hullStart;
                        do
                        {
                            if (_hullTri[e] == bl)
                            {
                                _hullTri[e] = a;
                                break;
                            }
                            e = _hullPrev[e];
                        } while (e != _hullStart);
                    }

                    Link(a, hbl);
                    Link(b, _halfedges[ar]);
                    Link(ar, bl);

                    var br = b0 + (b + 1) % 3;

                    // don't worry about hitting the cap: it can only happen on extremely degenerate input
                    if (i < _edgeStack.Length)
                    {
                        _edgeStack[i++] = br;
                    }
                }
                else
                {
                    if (i == 0) break;
                    a = _edgeStack[--i];
                }
            }

            return ar;
        }

        public Delaunator(double[] coords)
        {
            _coords = coords;
            _n = _coords.Length / 2;

            // arrays that will store the triangulation graph
            _maxTriangles = Math.Max(2 * _n - 5, 0);
            _triangles = new int[_maxTriangles * 3];
            _halfedges = new int[_maxTriangles * 3];

            // temporary arrays for tracking the edges of the advancing convex hull
            _hashSize = (int)Math.Ceiling(Math.Sqrt(_n));
            _hullPrev = new int[_n];
            _hullNext = new int[_n];
            _hullTri = new int[_n];
            _hullHash = new int[_hashSize];
            Array.Fill(_hullHash, -1);

            // temporary arrays for sorting points
            _ids = new int[_n];
            _dists = new double[_n];
        }

        public TriangulationInfo Build()
        {
            var minX = double.PositiveInfinity;
            var minY = double.PositiveInfinity;
            var maxX = double.NegativeInfinity;
            var maxY = double.NegativeInfinity;

            for (int i = 0; i < _n; i++)
            {
                var x = _coords[2 * i];
                var y = _coords[2 * i + 1];
                if (x < minX) minX = x;
                if (y < minY) minY = y;
                if (x > maxX) maxX = x;
                if (y > maxY) maxY = y;
                _ids[i] = i;
            }

            var cx = (minX + maxX) / 2.0;
            var cy = (minY + maxY) / 2.0;

            var minDist = double.PositiveInfinity;
            int i0 = 0;
            int i1 = 0;
            int i2 = 0;

            // pick a seed point close to the center
            for (int i = 0; i < _n; i++)
            {
                var d = DelaunatorHelpers.Dist(cx, cy, _coords[2 * i], _coords[2 * i + 1]);
                if (d < minDist)
                {
                    i0 = i;
                    minDist = d;
                }
            }

            var i0x = _coords[2 * i0];
            var i0y = _coords[2 * i0 + 1];

            minDist = double.PositiveInfinity;

            // find the point closest to the seed
            for (int i = 0; i < _n; i++)
            {
                if (i == i0) continue;
                var d = DelaunatorHelpers.Dist(i0x, i0y, _coords[2 * i], _coords[2 * i + 1]);
                if (d < minDist && d > 0)
                {
                    i1 = i;
                    minDist = d;
                }
            }

            var i1x = _coords[2 * i1];
            var i1y = _coords[2 * i1 + 1];

            var minRadius = double.PositiveInfinity;

            // find the third point which forms the smallest circumcircle with the first two
            for (int i = 0; i < _n; i++)
            {
                if (i == i0 || i == i1) continue;
                var r = DelaunatorHelpers.Circumradius(i0x, i0y, i1x, i1y, _coords[2 * i], _coords[2 * i + 1]);
                if (r < minRadius)
                {
                    i2 = i;
                    minRadius = r;
                }
            }

            var i2x = _coords[2 * i2];
            var i2y = _coords[2 * i2 + 1];

            if (minRadius == double.PositiveInfinity)
            {
                // order collinear points by dx (or dy if all x are identical)
                // and return the list as a hull
                for (int i = 0; i < _n; i++)
                {
                    var dx = (_coords[2 * i] - _coords[0]);
                    _dists[i] = dx > 0 ? dx : (_coords[2 * i + 1] - _coords[1]);
                }

                DelaunatorHelpers.Quicksort(_ids, _dists, 0, _n - 1);
                var hull = new int[_n];
                int j = 0;
                double d0 = double.NegativeInfinity;
                for (int i = 0; i < _n; i++)
                {
                    var id = _ids[i];
                    if (_dists[id] > d0)
                    {
                        hull[j++] = id;
                        d0 = _dists[id];
                    }
                }

                Array.Resize(ref hull, j);
                Array.Resize(ref _triangles, 0);
                Array.Resize(ref _halfedges, 0);
                return new TriangulationInfo(hull, _triangles, _coords, _halfedges);
            }

            // swap the order of the seed points for counter-clockwise orientation
            if (Orient2D.Robust(new[] { i0x, i0y }, new[] { i1x, i1y }, new[] { i2x, i2y }) > 0)
            {
                var i = i1;
                var x = i1x;
                var y = i1y;
                i1 = i2;
                i1x = i2x;
                i1y = i2y;
                i2 = i;
                i2x = x;
                i2y = y;
            }

            var center = DelaunatorHelpers.Circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
            _cx = center[0];
            _cy = center[1];

            for (int i = 0; i < _n; i++)
            {
                _dists[i] = DelaunatorHelpers.Dist(_coords[2 * i], _coords[2 * i + 1], _cx, _cy);
            }

            // sort the points by distance from the seed triangle circumcenter
            DelaunatorHelpers.Quicksort(_ids, _dists, 0, _n - 1);

            // set up the seed triangle as the starting hull
            _hullStart = i0;
            var hullSize = 3;

            _hullNext[i0] = _hullPrev[i2] = i1;
            _hullNext[i1] = _hullPrev[i0] = i2;
            _hullNext[i2] = _hullPrev[i1] = i0;

            _hullTri[i0] = 0;
            _hullTri[i1] = 1;
            _hullTri[i2] = 2;

            Array.Fill(_hullHash, -1);
            _hullHash[HashKey(i0x, i0y)] = i0;
            _hullHash[HashKey(i1x, i1y)] = i1;
            _hullHash[HashKey(i2x, i2y)] = i2;

            _trianglesLen = 0;
            AddTriangle(i0, i1, i2, -1, -1, -1);

            double xp = 0;
            double yp = 0;

            for (int k = 0; k < _ids.Length; k++)
            {
                var i = _ids[k];
                var x = _coords[2 * i];
                var y = _coords[2 * i + 1];

                // skip near-duplicate points
                if (k > 0 && Math.Abs(x - xp) <= _epsilon && Math.Abs(y - yp) <= _epsilon) continue;
                xp = x;
                yp = y;

                // skip seed triangle points
                if (i == i0 || i == i1 || i == i2) continue;

                // find a visible edge on the convex hull using edge hash
                int start = 0;
                var key = HashKey(x, y);
                for (int j = 0; j < _hashSize; j++)
                {
                    start = _hullHash[(key + j) % _hashSize];
                    if (start != -1 && start != _hullNext[start]) break;
                }

                start = _hullPrev[start];
                int e = start, q;
                q = _hullNext[e];
                while (Orient2D.Robust(
                    new[] { x, y },
                    new[] { _coords[2 * e], _coords[2 * e + 1] },
                    new[] { _coords[2 * q], _coords[2 * q + 1] }) <= 0)
                {
                    e = q;
                    if (e == start)
                    {
                        e = -1;
                        break;
                    }

                    q = _hullNext[e];
                }
                if (e == -1) continue; // likely a near-duplicate point; skip it

                // add the first triangle from the point
                var t = AddTriangle(e, i, _hullNext[e], -1, -1, _hullTri[e]);

                // recursively flip triangles from the point until they satisfy the Delaunay condition
                _hullTri[i] = Legalize(t + 2);
                _hullTri[e] = t; // keep track of boundary triangles on the hull
                hullSize++;

                // walk forward through the hull, adding more triangles and flipping recursively
                var n = _hullNext[e];
                q = _hullNext[n];
                while (Orient2D.Robust(new[] { x, y },
                    new[] { _coords[2 * n], _coords[2 * n + 1] },
                    new[] { _coords[2 * q], _coords[2 * q + 1] }) > 0)
                {
                    t = AddTriangle(n, i, q, _hullTri[i], -1, _hullTri[n]);
                    _hullTri[i] = Legalize(t + 2);
                    _hullNext[n] = n; // mark as removed
                    hullSize--;
                    n = q;
                    q = _hullNext[n];
                }

                // walk backward from the other side, adding more triangles and flipping
                if (e == start)
                {
                    q = _hullPrev[e];
                    while (Orient2D.Robust(
                        new[] { x, y },
                        new[] { _coords[2 * q], _coords[2 * q + 1] },
                        new[] { _coords[2 * e], _coords[2 * e + 1] }) > 0)
                    {
                        t = AddTriangle(q, i, e, -1, _hullTri[e], _hullTri[q]);
                        Legalize(t + 2);
                        _hullTri[q] = t;
                        _hullNext[e] = e; // mark as removed
                        hullSize--;
                        e = q;
                        q = _hullPrev[e];
                    }
                }

                // update the hull indices
                _hullStart = _hullPrev[i] = e;
                _hullNext[e] = _hullPrev[n] = i;
                _hullNext[i] = n;

                // save the two new edges in the hash table
                _hullHash[HashKey(x, y)] = i;
                _hullHash[HashKey(_coords[2 * e], _coords[2 * e + 1])] = e;
            }

            var finalHull = new int[hullSize];
            var ef = _hullStart;
            for (int i = 0; i < hullSize; i++)
            {
                finalHull[i] = ef;
                ef = _hullNext[ef];
            }

            Array.Resize(ref _triangles, _trianglesLen);
            Array.Resize(ref _halfedges, _trianglesLen);

            return new TriangulationInfo(finalHull, _triangles, _coords, _halfedges);
        }
    }
}
