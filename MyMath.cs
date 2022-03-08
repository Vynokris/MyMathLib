using System;
using System.Numerics;
using static System.Math;
using System.Diagnostics;

namespace MyMathLib
{
    // ---------- Colors ---------- //

    public class HSV
    {
        public float H { get; set; }
        public float S { get; set; }
        public float V { get; set; }

        public HSV() { H = 0; S = 0; V = 0; }
        public HSV(float h, float s, float v) { H = h; S = s; V = v; }

        // Convert an HSV color to RGB.
        public Color ToRGB(float alpha)
        {
            Color color = new Color( 0, 0, 0, alpha );
            float k = 0f, t = 0f;

            // Red channel
            k = (H + 5) % 6; t = 4 - k;
            k = (t < k) ? t : k;
            k = (k < 1) ? k : 1;
            k = (k > 0) ? k : 0;
            color.R = V - V * S * k;

            // Green channel
            k = (H + 3) % 6;  t = 4 - k;
            k = (t < k) ? t : k;
            k = (k < 1) ? k : 1;
            k = (k > 0) ? k : 0;
            color.G = V - V * S * k;

            // Blue channel
            k = (H + 1) % 6;  t = 4 - k;
            k = (t < k) ? t : k;
            k = (k < 1) ? k : 1;
            k = (k > 0) ? k : 0;
            color.B = V - V * S * k;

            return color;
        }
    }

    public class Color
    {
        public float R { get; set; }
        public float G { get; set; }
        public float B { get; set; }
        public float A { get; set; }

        public Color()                                   { R = 0; G = 0; B = 0; A = 1; }
        public Color(float r, float g, float b, float a) { R = r; G = g; B = b; A = a; }

        // Linear interpolation between two given colors.
        public void Lerp(float val, Color start, Color end)
        {
            R = start.R + val * (end.R - start.R);
            G = start.G + val * (end.G - start.G);
            B = start.B + val * (end.B - start.B);
            A = start.A + val * (end.A - start.A);
        }

        // Convert an RGB color (0 <= rgba <= 1) to HSV.
        public HSV ToHSV()
        {
            HSV hsv = new HSV();

            float minV = Min(Min(R, G), B);
            float maxV = Max(Max(R, G), B);
            float diff = maxV - minV;

            float r = R, g = G, b = B;

            // Set Value.
            hsv.V = maxV;

            // If max and min are the same, return.
            if (diff< 0.00001f) return new HSV(0, 0, hsv.V);

            // Set Saturation.
            if (maxV > 0) hsv.S = diff / maxV;
            else          return new HSV(0, 0, hsv.V);

            // Set Hue.
            if      (r >= maxV) hsv.H = (g - b) / diff;
            else if (g >= maxV) hsv.H = 2.0f + (b - r) / diff;
            else if (b >= maxV) hsv.H = 4.0f + (r - g) / diff;

            // Keep Hue above 0.
            if (hsv.H < 0) hsv.H += 2 * (float)PI;

            return hsv;
        }

        // Shifts the hue of the given color.
        public Color Shift(float hue)
        {
            HSV hsv = this.ToHSV();
            hsv.H += hue;
            if (hsv.H >= 2 * (float)PI) hsv.H -= 2 * (float)PI;
            else if (hsv.H < 0)         hsv.H += 2 * (float)PI;
            Color result = hsv.ToRGB(A);

            return result;
        }
    }



    // ---------- Arithmetic ---------- //

    public static class Arithmetic
    {
        // Rounds the given value to the nearest int.
        public static int RoundInt(float val)   { return (int)Round(val);  }

        // Rounds down the given value.
        public static int FloorInt(float val)   { return (int)Floor(val);  }

        // Rounds up the given value.
        public static int CeilInt(float val)    { return (int)Ceiling(val); }

        // Returns the sqare power of the given value.
        public static float SqPow(float val)    { return val * val;  }

        // Returns 1 if the given value is positive or null, and -1 if it is negative.
        public static int SignOf(float val)     { if (val == 0) return 1; return (int)val / Abs((int)val); }

        // Converts the given angle from degrees to radians.
        public static float Deg2Rad(float val)  { return val * ((float)PI / 180f); }

        // Converts the given angle from radians to degrees.
        public static float RadToDeg(float rad) { return rad * (180f / (float)PI); }

        // Clamps the given value to be inferior or equal to the maximum value.
        public static float ClampUnder(float val, float max)        { if (val > max) val = max; return val; }

        // Clamps the given value to be superior or equal to the minimum value.
        public static float ClampAbove(float val, float min)        { if (val < min) val = min; return val; }

        // Compute linear interpolation between start and end for the parameter val (if 0 <= val <= 1: start <= return <= end).
        public static float Lerp(float val, float start, float end) { return start + val* (end - start); }

        // Compute the linear interpolation factor that returns val when lerping between start and end.
        public static float GetLerp(float val, float start, float end)
        {
            if (end - start != 0) return (val - start) / (end - start);
            return 0;
        }

        // Remaps the given value from one range to another.
        public static float Remap(float val, float inputStart, float inputEnd, float outputStart, float outputEnd)
        {
            return outputStart + (val - inputStart) * (outputEnd - outputStart) / (inputEnd - inputStart);
        }

        // Returns true if the given number is a power of 2.
        public static bool IsPowerOf2(int val)      { return val == (int)Pow(2, (int)(Log(val) / Log(2))); }

        // Returns the closest power of 2 that is inferior or equal to val.
        public static int GetPowerOf2Under(int val) {  return (int)Pow(2, (int)(Log(val) / Log(2))); }

        // Returns the closest power of 2 that is superior or equal to val.
        public static int GetPowerOf2Above(int val)
        {
            if (IsPowerOf2(val)) return (int)Pow(2, (int)(Log(val) / Log(2)));
            else                 return (int)Pow(2, (int)(Log(val) / Log(2)) + 1);
        }

        // Blend between two HSV colors.
        public static HSV BlendHSV(HSV color0, HSV color1)
        {
            Vector2 totalVec = Geometry2D.Vector2FromAngle(color0.H, 1)
                             + Geometry2D.Vector2FromAngle(color1.H, 1);

            float avgHue = totalVec.GetAngle();
            float avgSat = (color0.S + color1.S) / 2;
            float avgVal = (color0.V + color1.V) / 2;

            return new HSV(avgHue, avgSat, avgVal);
        }
    }



    // ---------- Geometry2D ---------- //

    public static class Geometry2D
    {
        // ---- Vector2 ---- //

        // Returns a null vector.
        public static Vector2 Vector2Zero()                             { return new Vector2(0, 0); }

        // Creates a vector from two values.
        public static Vector2 Vector2Create(float X, float Y)           { return new Vector2(X, Y); }

        // Creates a vector from one point to another.
        public static Vector2 Vector2FromPoints(Vector2 p0, Vector2 p1) { return new Vector2(p1.X - p0.X, p1.Y - p0.Y); }

        // Creates a vector given an angle and a length.
        public static Vector2 Vector2FromAngle(float rad, float length) { return new Vector2((float)Cos(rad) * length,
                                                                                             (float)Sin(rad) * length); }
        // Creates a vector from a segment.
        public static Vector2 Vector2FromSegment(Segment2 seg)          { return Vector2FromPoints(seg.A, seg.B); }
        
        // Vector dot product.
        public static float Dot(this Vector2 v0, Vector2 v1) { return (v0.X * v1.X) + (v0.Y * v1.Y); }

        // Vector cross product.
        public static float Cross(this Vector2 v0, Vector2 v1) { return (v0.X * v1.Y) - (v0.Y * v1.X); }

        // Returns a normalized copy of the vector.
        public static Vector2 GetNormalized(this Vector2 v)             { return v / v.Length(); }

        // Normalizes the vector so that its length is 1.
        public static void Normalize(this ref Vector2 v)                { v = v.GetNormalized(); }

        // Modifies the length of the given vector to correspond to the given value.
        public static void SetLength(this ref Vector2 v, float length)  { v = v.GetNormalized() * length; }

        // Returns a vector width the destination's size and the source's signs.
        public static Vector2 GetCopiedSign(this Vector2 dest, Vector2 source)
        {
            return new Vector2((float)Math.CopySign((float)dest.X, (float)source.X),
                               (float)Math.CopySign((float)dest.Y, (float)source.Y));
        }

        // Copies the signs from the source vector to the destination vector.
        public static void CopySign(this ref Vector2 dest, Vector2 source)
        {
            dest = dest.GetCopiedSign(source);
        }

        // Returns the normal of the vector.
        public static Vector2 GetNormal(this Vector2 v) { return new Vector2(-v.Y, v.X); }

        // Interprets the vector as a point and returns the distance to another point.
        public static float GetDistanceFromPoint(this Vector2 p0, Vector2 p1) { return (p1 - p0).Length(); }

        // Returns the angle (in radians) of the vector.
        public static float GetAngle(this Vector2 v)
        {
            return (float)Math.CopySign((float)Acos((double)v.GetNormalized().X),
                                        (float)Asin((double)v.GetNormalized().Y));
        }

        // Returns the angle (in radians) between two vectors.
        public static float GetAngleWithVector(this Vector2 v0, Vector2 v1)
        {
            float v0angle = v0.GetAngle();
            float v1angle = v1.GetAngle();
            return (v0angle >= v1angle ? (v0angle - v1angle) : (v1angle - v0angle));
        }

        // Rotates the given vector by the given angle.
        public static void Rotate(this ref Vector2 v, float angle)
        {
            float vLength = v.Length();
            float vAngle = v.GetAngle();
            v = new Vector2((float)Cos(vAngle + angle) * vLength, (float)Sin(vAngle + angle) * vLength);
        }
        
        // Calculates linear interpolation for a value from a start point to an end point.
        public static Vector2 Point2Lerp(float val, Vector2 start, Vector2 end)
        {
            return new Vector2(Arithmetic.Lerp(val, start.X, end.X),
                               Arithmetic.Lerp(val, start.Y, end.Y));
        }


        // ---- Segment2 ---- //
        public class Segment2
        {
            public Vector2 A, B;

            // Null Segment.
            public Segment2() { A = Vector2Zero(); B = Vector2Zero(); }

            // Segment from points.
            public Segment2(Vector2 a, Vector2 b) { A = a; B = b; }

            // Segment from point and vector.
            public Segment2(Vector2 origin, Vector2 vec, bool vector) { A = origin; B = origin + vec; }

            // Returns the center of mass of the Segment2.
            public Vector2 GetCenterOfMass() { return Vector2Create((A.X + B.X) / 2, (A.Y + B.Y) / 2); }

            // Returns the number of sides of the Segment2.
            public int GetSidesNum() { return 1; }

            // Returns the side of the Segment2 that corresponds to the given index.
            public Segment2 GetSide(int index)
            {
                Debug.Assert(0 <= index && index < 1);
                return this;
            }

            // Returns the number of vertices of the Segment2.
            public int GetVerticesNum() { return 2; }

            // Returns the vertex of the Segment2 that corresponds to the given index.
            public Vector2 GetVertex(int index)
            {
                Debug.Assert(0 <= index && index< 2);

                return index switch {
                    0 => A, 
                    1 => B, 
                    _ => Vector2Zero(),
                };
            }

            // Moves the Segment2 by the given vector.
            public void Move(Vector2 vec) { A += vec; B += vec; }
        }


        // ---- Triangle2 ---- //
        public class Triangle2
        {
            public Vector2 A, B, C;

            // Null triangle.
            public Triangle2() { A = Vector2Zero(); B = Vector2Zero(); C = Vector2Zero(); }

            // Triangle from points.
            public Triangle2(Vector2 a, Vector2 b, Vector2 c) { A = a; B = b; C = c; }

            // Returns the center of mass of the triangle.
            public Vector2 GetCenterOfMass() { return Vector2Create((A.X + B.X + C.X) / 3,
                                                                    (A.Y + B.Y + C.Y) / 3); }

            // Returns the number of sides of the triangle.
            public int GetSidesNum() { return 3; }

            // Returns the side of the triangle that corresponds to the given index.
            public Segment2 GetSide(int index)
            {
                Debug.Assert(0<= index && index< 3);

                return index switch
                {
                    0 => new Segment2(A, B),
                    1 => new Segment2(B, C),
                    2 => new Segment2(C, A),
                    _ => new Segment2(Vector2Zero(), Vector2Zero()),
                };
            }

            // Returns the number of vertices of the triangle.
            public int GetVerticesNum() { return 3; }

            // Returns the vertex of the triangle that corresponds to the given index.
            public Vector2 GetVertex(int index)
            {
                Debug.Assert(0 <= index && index< 3);

                return index switch
                {
                    0 => A,
                    1 => B,
                    2 => C,
                    _ => Vector2Zero(),
                };
            }

            // Moves the triangle by the given vector.
            public void Move(Vector2 vec) { A += vec; B += vec; C += vec; }
        }


        // ---- Rectangle2 ---- //
        public class Rectangle2
        {
            public Vector2 O;
            public float W, H;

            // Null rectangle.
            public Rectangle2() { O = Vector2Zero(); W = 0; H = 0; }

            // Rectangle from posX posY, width and height.
            public Rectangle2(float x, float y, float w, float h) { O = Vector2Create(x, y); W = w; H = h; }

            // Returns the center of mass of the rectangle.
            public Vector2 GetCenterOfMass() { return Vector2Create(O.X + W / 2, O.Y + H / 2); }

            // Returns the number of sides of the rectangle.
            public int GetSidesNum() { return 4; }

            // Returns the side of the rectangle that corresponds to the given index.
            public Segment2 GetSide(int index)
            {
                Debug.Assert(0 <= index && index < 4);

                return index switch
                {
                    0 => new Segment2(Vector2Create(O.X + W, O.Y), O),
                    1 => new Segment2(O, Vector2Create(O.X, O.Y + H)),
                    2 => new Segment2(Vector2Create(O.X, O.Y + H), Vector2Create(O.X + W, O.Y + H)),
                    3 => new Segment2(Vector2Create(O.X + W, O.Y + H), Vector2Create(O.X + W, O.Y)),
                    _ => new Segment2(Vector2Zero(), Vector2Zero()),
                };
            }

            // Returns the number of vertices of the rectangle.
            public int GetVerticesNum() { return 4; }

            // Returns the vertex of the rectangle that corresponds to the given index.
            public Vector2 GetVertex(int index)
            {
                Debug.Assert(0 <= index && index< 4);

                return index switch
                {
                    0 => Vector2Create(O.X + W, O.Y),
                    1 => O,
                    2 => Vector2Create(O.X, O.Y + H),
                    3 => Vector2Create(O.X + W, O.Y + H),
                    _ => Vector2Zero(),
                };
            }

            // Moves the rectangle by the given vector.
            public void Move(Vector2 vec) { O += vec; }
        }


        // ---- Polygon2 ---- //
        public class Polygon2
        {
            public Vector2 O;
            public float Radius, Rot;
            public int Sides;

            // Null polygon.
            public Polygon2() { O = Vector2Zero(); Radius = 0; Rot = 0; Sides = 3; }

            // Polygon from origin, radius, rotation and number of sides.
            public Polygon2(float x, float y, float radius, float rotation, int sides)
            {
                O = Vector2Create(x, y); Radius = radius; Rot = rotation; Sides = sides;
            }

            // Returns the center of mass of the polygon.
            public Vector2 GetCenterOfMass() { return O; }

            // Returns the number of sides of the polygon.
            public int GetSidesNum() { return Sides; }

            // Returns the side of the polygon that corresponds to the given index.
            public Segment2 GetSide(int index)
            {
                Debug.Assert(0 <= index && index < Sides);

                float corner_angle = Arithmetic.Deg2Rad(360 / Sides);
                float angle_offset = (float)PI / 2 + (index * corner_angle);
                Vector2 poly_point_a = O + Vector2FromAngle(angle_offset + Rot, Radius);
                Vector2 poly_point_b = O + Vector2FromAngle(angle_offset + corner_angle + Rot, Radius);

                return new Segment2(poly_point_a, poly_point_b);
            }

            // Returns the number of vertices of the polygon.
            public int GetVerticesNum() { return Sides; }

            // Returns the vertex of the polygon that corresponds to the given index.
            public Vector2 GetVertex(int index)
            {
                Debug.Assert(0 <= index && index < Sides);

                float corner_angle = Arithmetic.Deg2Rad(360 / Sides);
                float angle_offset = (float)PI / 2 + (index * corner_angle);
                return O + Vector2FromAngle(angle_offset + Rot, Radius);
            }

            // Moves the polygon by the given vector.
            public void Move(Vector2 vec) { O += vec; }
        }


        // ---- Circle2 ---- //
        public class Circle2
        {
            public Vector2 O;
            float Radius;

            // Null circle.
            public Circle2() { O = Vector2Zero(); Radius = 0; }

            // Circle from position and radius.
            public Circle2(float x, float y, float radius) { O = Vector2Create(x, y); Radius = radius; }

            // Returns the center of mass of the circle.
            public Vector2 GetCenterOfMass() { return O; }

            // Returns the number of sides of the circle.
            public int GetSidesNum() { return 1; }

            // Does nothing and returns a null Segment2.
            public Segment2 GetSide(int index) { return new Segment2(); }

            // Returns the number of vertices of the circle.
            public int GetVerticesNum() { return 0; }

            // Does nothing and returns a null vector.
            public Vector2 GetVertex(int index) { return Vector2Zero(); }

            // Moves the circle by the given vector.
            public void Move(Vector2 vec) { O += vec; }
        }
    }



    // ---------- Geometry3D ---------- //

    public struct Vertex
    {
        Vector3 pos;
        Vector3 normal;
        Color   color;
        Vector2 uv;
    }

    public static class Geometry3D
    {
        // Returns the coordinates of a point on a sphere of radius r, using the given angles.
        public static Vector3 GetSphericalCoords(float r, float theta, float phi)
        {
            return new Vector3(r * (float)Sin(theta) * (float)Cos(phi),
                               r * (float)Cos(theta),
                               r * (float)Sin(theta) * (float)Sin(phi));
        }


        // ---- Vector3 ---- //

        // Null vector.
        public static Vector3 Vector3Zero() { return new Vector3(0, 0, 0); }

        // Vector from 3 coordinates.
        public static Vector3 Vector3Create(float x, float y, float z) { return new Vector3(x, y, z); }

        // Vector from points.
        public static Vector3 Vector3FromPoints(Vector3 p0, Vector3 p1) { return new Vector3(p1.X - p0.X, p1.Y - p0.Y, p1.Z - p0.Z); }

        // Vector from angle.
        public static Vector3 Vector3FromAngle(float theta, float phi, float length)
        {
            return new Vector3(length * (float)Sin(theta) * (float)Cos(phi),
                               length * (float)Cos(theta),
                               length * (float)Sin(theta) * (float)Sin(phi));
        }

        // Vector from segment.
        // TODO.

        // Vector dot product.
        public static float Dot(this Vector3 v0, Vector3 v1) { return (v0.X * v1.X) + (v0.Y * v1.Z) + (v0.Z * v1.Z); }

        // Vector cross product.
        public static Vector3 Cross(this Vector3 v0, Vector3 v1) { return new Vector3((v0.Y * v1.Z - v0.Z * v1.Y), (v0.Z * v1.X - v0.X * v1.Z), (v0.X * v1.Y - v0.Y * v1.X)); }

        // Normalizes the vector so that its length is 1.
        public static void Normalize(this ref Vector3 v) { v /= v.Length(); }

        // Returns a normalized copy of the vector.
        public static Vector3 GetNormalized(this Vector3 v) { return v / v.Length(); }

        // Sets the length of the vector to the given value.
        public static void SetLength(this ref Vector3 v, float length) { v = v.GetNormalized() * length; }

        // Sets the destination's sign to those of the source.
        public static void CopySign(this ref Vector3 dest, Vector3 source) { dest = dest.GetCopiedSign(source); }

        // Returns a new vector width the destination's size and the source's signs.
        public static Vector3 GetCopiedSign(this Vector3 dest, Vector3 source)
        {
            return new Vector3((float)Math.CopySign(dest.X, source.X),
                               (float)Math.CopySign(dest.Y, source.Y),
                               (float)Math.CopySign(dest.Z, source.Z));
        }

        // Interprets the vector as a point and returns the distance to another point.
        public static float GetDistanceFromPoint(this Vector3 p0, Vector3 p1) { return Vector3FromPoints(p0, p1).Length(); }

        // Returns the angle (in radians) of the given vector.
        public static float GetAngleTheta(this Vector3 v)
        {
            return (2*(float)PI - Geometry2D.Vector2Create(v.X, v.Z).GetAngle() + (float)PI/2) % 2*(float)PI;
        }
        public static float GetAnglePhi  (this Vector3 v) 
        {
            return (2 * (float)PI - Geometry2D.Vector2Create(Geometry2D.Vector2Create(v.X, v.Z).Length(), v.Y).GetAngle()) - 2 * (float)PI;
        }

        // Returns the angle (in radians) between two vectors.
        public static float GetAngleThetaWithVec3(this Vector3 v0, Vector3 v1)
        {
            float v0angle = v0.GetAngleTheta();
            float v1angle = v1.GetAngleTheta();
            return (v0angle >= v1angle ? (v0angle - v1angle) : (v1angle - v0angle));
        }
        public static float GetAnglePhiWithVec3(this Vector3 v0, Vector3 v1)
        {
            float v0angle = v0.GetAnglePhi();
            float v1angle = v1.GetAnglePhi();
            return (v0angle >= v1angle? (v0angle - v1angle) : (v1angle - v0angle));
        }

        // Rotates the given vector by the given angle.
        public static void Rotate(this ref Vector3 v, float theta, float phi) 
        {
            v = Vector3FromAngle(v.GetAngleTheta() + theta, v.GetAnglePhi() + phi, v.Length());
        }

        // Creates a Vector4 from this vector.
        public static Vector4 ToVec4(this Vector3 v) { return Vector4Create(v.X, v.Y, v.Z, 1); }

        // Calculates linear interpolation for a value from a start point to an end point.
        public static Vector3 Point3Lerp(float val, Vector3 start, Vector3 end)
        {
        return new Vector3(Arithmetic.Lerp(val, start.X, end.X),
                           Arithmetic.Lerp(val, start.Y, end.Y),
                           Arithmetic.Lerp(val, start.Z, end.Z));
        }


        // ---- Vector4 ---- //

        // Null vector.
        public static Vector4 Vecto4Zero() { return new Vector4(0, 0, 0, 1); }

        // Vector from coordinates.
        public static Vector4 Vector4Create(float x, float y, float z, float w) { return new Vector4(x, y, z, w); }

        // Vector from 2 points.
        public static Vector4 Vector4FromPoints(Vector4 p0, Vector4 p1, float w) { return new Vector4(p1.X - p0.X, p1.Y - p0.Y, p1.Z - p0.Z, w); }

        // Vector4 from vector3.
        public static Vector4 Vector4FromVec3(Vector3 v, float w) { return new Vector4(v.X, v.Y, v.Z, w); }

        // Vector from segment.
        // TODO.

        // Vector from angle.
        public static Vector4 Vector4FromAngle(float theta, float phi, float length, float w)
        {
            return new Vector4(length * (float)Sin(theta) * (float)Cos(phi),
                               length * (float)Cos(theta),
                               length * (float)Sin(theta) * (float)Sin(phi),
                               w);
        }

        // Vector dot product.
        public static float Dot(this Vector4 v0, Vector4 v1) { return v0.ToVec3(false).Dot(v1.ToVec3(false)); }

        // Vector cross product.
        public static Vector3 Cross(this Vector4 v0, Vector4 v1) { return v0.ToVec3(false).Cross(v1.ToVec3(false)); }

        // Homogenizes the vector4 by dividing it by w.
        public static Vector4 Homogenize(this Vector4 v)
        {
            return new Vector4(v.X / v.W, v.Y / v.W, v.Z / v.W, 1f);
        }
        
        // Normalizes the given vector so that its length is 1.
        public static Vector4 Normalize(this Vector4 v) 
        {
            return new Vector4(v.X / v.Length(), v.Y / v.Length(),
                               v.Z / v.Length(), v.W / v.Length());
        }

        // Negates both of the coordinates of the given vector.
        public static Vector4 Negate(this Vector4 v) { return new Vector4(-v.X, -v.Y, -v.Z, v.W); }

        // Copies the signs from the source vector to the destination vector.
        public static Vector4 Copysign(this Vector4 dest, Vector4 src) 
        {
            return new Vector4((float)Math.CopySign(dest.X, src.X), (float)Math.CopySign(dest.Y, src.Y),
                               (float)Math.CopySign(dest.Z, src.Z), (float)Math.CopySign(dest.W, src.W));
        }

        // Interprets the vector as a point and returns the distance to another point.
        public static float GetDistanceFromPoint(this Vector4 p0, Vector4 p1)
        {
            return new Vector4(p1.X - p0.X, p1.Y - p0.Y, p1.Z - p0.Z, p0.W).Length();
        }

        // Returns the angle (in radians) of the given vector.
        public static float GetAngleTheta(this Vector4 v) { return (float)Acos(v.Z / v.Length()); }

        public static float GetAnglePhi  (this Vector4 v)
        { 
            if (v.X < 0) return (float)Atan(v.Y / v.X);
            if (v.X > 0) return (float)Atan(v.Y / v.X) + (float)PI;
            return (float)PI / 2;
        }

        // Returns the angle (in radians) between two vectors.
        public static float GetAngleThetaWithVec4(this Vector4 v1, Vector4 v2)
        {
            float v1Angle = v1.GetAngleTheta();
            float v2Angle = v2.GetAngleTheta();
            return (v1Angle >= v2Angle ? (v1Angle - v2Angle) : (v2Angle - v1Angle));
        }

        public static float GetAnglePhiWithVec4(this Vector4 v1, Vector4 v2)
        {
            float v1Angle = v1.GetAnglePhi();
            float v2Angle = v2.GetAnglePhi();
            return (v1Angle >= v2Angle ? (v1Angle - v2Angle) : (v2Angle - v1Angle));
        }

        // Rotates the given vector by the given angle.
        public static void Rotate(this ref Vector4 v, float theta, float phi)
        {
            v = Vector4FromAngle(theta,  phi, v.Length(), v.W);
        }

        // Creates a Vector3 from this vector.
        public static Vector3 ToVec3(this Vector4 v, bool homogenizeVec)
        {
            if (homogenizeVec)
                return Vector3Create(v.X/v.W, v.Y/v.W, v.Z/v.W);
            else
                return Vector3Create(v.X, v.Y, v.Z);
        }
    }
}
