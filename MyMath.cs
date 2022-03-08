using System;

namespace MyMathLib
{
    public static class Arithmetic
    {
        public static int RoundInt(ref float val)  { return (int)Math.Round(val);  }
        public static int FloorInt(ref float val)  { return (int)Math.Floor(val);  }
        public static int CeilInt(ref float val)   { return (int)Math.Ceiling(val); }
        public static float SqPow(ref float val)   { return val * val;  }
        public static int SignOf(ref float val)    { if (val == 0) return 1; return (int)val / Math.Abs((int)val); }
        public static float Deg2Rad(ref float val) { return val * ((float)Math.PI / 180f); }
    }
}
