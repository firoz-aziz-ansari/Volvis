/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class VoxelGradient {

    public float x, y, z;
    public float mag;
    
    public VoxelGradient() {
        x = y = z = mag = 0.0f;
    }
    
    public VoxelGradient(float gx, float gy, float gz) {
        x = gx;
        y = gy;
        z = gz;
        mag = (float) Math.sqrt(x*x + y*y + z*z);
    }
    
    public float getMag()
    {
       return this.mag;
    }
    public double[] getNormal() {
        mag = (mag == 0.0) ? (float)0.0001 : mag;
        return new double[] {x/mag, y/mag, z/mag};
    }
}
