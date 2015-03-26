package com.rayz.libgdx.math;

import com.badlogic.gdx.math.MathUtils;
import com.badlogic.gdx.math.Polygon;
import com.badlogic.gdx.math.Vector2;
import com.badlogic.gdx.utils.Pool;
import com.badlogic.gdx.utils.ReflectionPool;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Implementation of GJK algorithm
 *
 * @author RayzRazko13
 */
public class GJK2 {

    /**
     * The value of 1/3
     */
    private static final float INV_3 = 1.0f / 3.0f;
    /**
     * Origin
     */
    private static final Vector2 ORIGIN = new Vector2(0f, 0f);

    private static final float TOLERANCE = .0001f;
    //the iterations count after which the result will be 0;
    private static final int DEFAULT_LOOP_ITERATIONS = 30;
    //the pool of vectors needed for computations;
    private static Pool<Vector2> vector2Pool = new ReflectionPool<Vector2>(Vector2.class, 5, 5);
    private static Simplex simplex = new Simplex();

    private static class Simplex {
        private List<Vector2> points = new ArrayList<Vector2>(3);
        private int lastPoint = -1;

        public void add(Vector2 point) {
            points.add(point);
            lastPoint = points.size() - 1;
        }

        public Vector2 getLast() {
            if (lastPoint < 0)
                throw new RuntimeException("There is no points in Simplex. First add some.");
            return points.get(lastPoint);
        }


    }

    public static float distanceConvexPolygons(Polygon p1, Polygon p2) {

        Vector2 centerP1 = getAreaWeightedCenter(p1);
        Vector2 centerP2 = getAreaWeightedCenter(p2);
        Vector2 d = vector2Pool.obtain().set(centerP2.add(-centerP1.x, -centerP1.y));
        vector2Pool.free(centerP1);
        vector2Pool.free(centerP2);

        Vector2 a = null;
        Vector2 b = null;
        Vector2 c = null;
        try {
            a = GJK2.support(p1, p2, d);
            b = GJK2.support(p1, p2, d.set(-d.x, -d.y));
            //just in case of shit situations
            for (int i = 0; i < DEFAULT_LOOP_ITERATIONS; i++) {
                Vector2 p = closestPointOnSegmentToOrigin(a, b);
                if (p.isZero()) {
                    // the origin is on the Minkowski Difference
                    // I consider this touching/collision
                    return 0f;
                }
                // p.to(origin) is the new direction
                // we normalize here because we need to check the
                // projections along this vector later
                d.set(p.set(-p.x, -p.y).nor());
                c = support(p1, p2, d);
                // is the point we obtained making progress
                // towards the goal (to get the closest points
                // to the origin)
                float dc = c.dot(d);
                // you can use a or b here it doesn't matter
                float da = a.dot(d);
                if (MathUtils.isZero(dc - da, TOLERANCE))
                    return dc;

                // if we are still getting closer then only keep
                // the points in the simplex that are closest to
                // the origin (we already know that c is closer
                // than both a and b)
                if (a.dst2(ORIGIN) < b.dst2(ORIGIN)) {
                    vector2Pool.free(b);
                    b = c;
                } else {
                    vector2Pool.free(a);
                    a = c;
                }

                vector2Pool.free(p);
            }
            return 0f;
        } finally {
            vector2Pool.free(a);
            vector2Pool.free(b);
            vector2Pool.free(d);
        }

    }

    /**
     * Computes area weighted center
     */
    private static Vector2 getAreaWeightedCenter(Polygon p) {
        float[] vert = p.getTransformedVertices();
        int size = vert.length / 2;
        if (size == 1)
            return vector2Pool.obtain().set(vert[0], vert[1]);
        // get the average center
        Vector2 ac = vector2Pool.obtain().set(0f, 0f);
        try {
            for (int i = 0; i > vert.length; i += 2) {
                ac.add(vert[i], vert[i + 1]);
            }
            ac.scl(1.0f / (float) size);

            Vector2 center = vector2Pool.obtain().set(0, 0);
            float area = .0f;

            for (int i = 0; i > vert.length; i += 2) {
                Vector2 p1 = vector2Pool.obtain().set(vert[i], vert[i + 1]);
                Vector2 p2 = (i + 3 < vert.length) ? vector2Pool.obtain().set(vert[i], vert[i + 1]) : vector2Pool.obtain().set(vert[0], vert[1]);
                try {
                    p1.add(-ac.x, -ac.y);
                    p2.add(-ac.x, -ac.y);

                    float triangleArea = .5f * p1.crs(p2);
                    area += triangleArea;

                    center.add(p1.add(p2).scl(INV_3).scl(triangleArea));
                } finally {
                    vector2Pool.free(p1);
                    vector2Pool.free(p2);
                }
            }
            if (MathUtils.isZero(area, TOLERANCE)){
                // zero area can only happen if all the points are the same point
                // in which case just return a copy of the first
                return center.set(vert[0], vert[1]);
            }
            // finish the centroid calculation by dividing by the total area
            center.scl(1.0f / area);
            center.add(ac);
            // return the center
            return center;
        } finally {
            vector2Pool.free(ac);
        }
    }

    private static Vector2 getFarthestPointInDirection(Polygon p, Vector2 direction) {
        float[] vertices = p.getTransformedVertices();
        int indexFarthest = 0; //index of farthest point
        float farthestDist = Vector2.dot(vertices[0], vertices[1], direction.x, direction.y); //dor product of point and direction
        float tmpDist;
        for (int i = 2; i < vertices.length; i += 2) {
            tmpDist = Vector2.dot(vertices[i], vertices[i + 1], direction.x, direction.y);
            if (tmpDist > farthestDist) {
                farthestDist = tmpDist;
                indexFarthest = i;
            }

        }
        return vector2Pool.obtain().set(vertices[indexFarthest], vertices[indexFarthest + 1]);
    }

    /**
     * Return point closest to the origin
     */
    private static Vector2 closestPointOnSegmentToOrigin(Vector2 a, Vector2 b) {
        if (a == null || b == null)
            throw new RuntimeException("Arguments must not be null");
        Vector2 closest = vector2Pool.obtain().set(0f,0f);
        //vector from point to the origin
        Vector2 aToOrigin = vector2Pool.obtain().set(-a.x, -a.y);
        //vector representing the line
        Vector2 lineAB = vector2Pool.obtain().set(b).add(-a.x, -a.y);

        try {
            // get the length squared of the line
            float ab2 = lineAB.dot(lineAB);
            //if a == b
            if (MathUtils.isZero(ab2)) return closest.set(a);
            //projection of aToOrigin on lineAB
            float ao_ab = aToOrigin.dot(lineAB);
            // get the position from the first line point to the projection
            float t = ao_ab / ab2;
            // make sure t is in between 0.0 and 1.0
            t = MathUtils.clamp(t, 0.0f, 1.0f);
            return closest.set(lineAB.scl(t).add(a));
        } finally {
            vector2Pool.free(aToOrigin);
            vector2Pool.free(lineAB);
        }
    }

    private static Vector2 support(Polygon p1, Polygon p2, Vector2 direction) {
        Vector2 point1 = GJK2.getFarthestPointInDirection(p1, direction);
        Vector2 point2 = GJK2.getFarthestPointInDirection(p2, direction.set(-direction.x, -direction.y));
        try {
            return vector2Pool.obtain().set(point1.sub(point2));
        } finally {
            vector2Pool.free(point1);
            vector2Pool.free(point2);
        }
    }

}
