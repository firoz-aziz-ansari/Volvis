/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    int opt =1;
    boolean shading = false;
    double[] viewVec;
    double[] uVec;
    double[] vVec;
    double[] pixelCoord;
    double[] volumeCenter;
    TFColor voxelColor;
    double max;
    int imageCenter;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    double getVoxel(double[] coord) {
        //Trilinear Interpolation https://gist.github.com/begla/1019993
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];
        
        // calculating floor and ceil
        int x0 = (int) Math.floor(x);
        int x1 = (int) Math.min(Math.ceil(x), volume.getDimX() - 1);
        int y0 = (int) Math.floor(y);
        int y1 = (int) Math.min(Math.ceil(y), volume.getDimY() - 1);
        int z0 = (int) Math.floor(z);
        int z1 = (int) Math.min(Math.ceil(z), volume.getDimZ() - 1);
        
         if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ() || x0 >= volume.getDimX() || x0 < 0 || y0 >= volume.getDimY() || y0 < 0 || z0 >= volume.getDimZ() || z0<0 || x1<0 || x1 > volume.getDimX() || y1 < 0 || y1 > volume.getDimY() || z1<0 || z1> volume.getDimZ()) {
            return 0;
        }
        
        // calculating the 8 corners of the cube
        double s000 = volume.getVoxel(x0, y0, z0);
        double s001 = volume.getVoxel(x0, y0, z1);
        double s010 = volume.getVoxel(x0, y1, z0);
        double s100 = volume.getVoxel(x1, y0, z0);
        double s101 = volume.getVoxel(x1, y0, z1);
        double s011 = volume.getVoxel(x0, y1, z1);
        double s110 = volume.getVoxel(x1, y1, z0);
        double s111 = volume.getVoxel(x1, y1, z1);
        
        // x axis
        double s00 = linearInterpolate(x, x0, x1, s000, s100);
        double s01 = linearInterpolate(x, x0, x1, s001, s101);
        double s10 = linearInterpolate(x, x0, x1, s010, s110);
        double s11 = linearInterpolate(x, x0, x1, s011, s111);
        
        // y axis
        double s0 = linearInterpolate(y, y0, y1, s00, s10);
        double s1 = linearInterpolate(y, y0, y1, s01, s11);
        
        // z axis
        double s = linearInterpolate(z, z0, z1, s0, s1);
        return s;
    }
    
     public static double linearInterpolate(double x, double x0, double x1, double v0, double v1){
           double alpha = (x-x0)/(x1-x0);
           return (1 - alpha)*v0 + alpha*v1;
    }


    void slicer(double[] viewMatrix) {
        
        
        int res =1;
        if(interactiveMode) res=3;
        
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        viewVec = new double[3];
        uVec = new double[3];
        vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        imageCenter = image.getWidth() / 2;

        pixelCoord = new double[3];
        volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        max = volume.getMaximum();
        voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j=j+res) {
            for (int i = 0; i < image.getWidth(); i=i+res) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

                double val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int l = 0; l <res; l++) {
                    for (int m= 0; m <res; m++) {
                        if (m + i < image.getHeight() && m + i >= 0 && l + j < image.getWidth() && l + j >= 0) {
                            image.setRGB(m + i, l + j,pixelColor);
                        }
                    }
                }
            }
        }

    }
    
    void mip(double[] viewMatrix){
         // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        int res =1;
        if(interactiveMode) res=3;
        
        viewVec = new double[3];
        uVec = new double[3];
        vVec = new double[3];     
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        imageCenter = image.getWidth() / 2;
        pixelCoord = new double[3];
        volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        max = volume.getMaximum();
        
        for (int j = 0; j < image.getHeight(); j+= res) {
            for (int i = 0; i < image.getWidth(); i+= res) {
                
                int castRay = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()); 
                double maxInt = 0;
                
                for (int k = 0; k < castRay - 1; k++) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) +(k - imageCenter) * viewVec[0] + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + (k - imageCenter) * viewVec[1] + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) +  (k - imageCenter)*viewVec[2] + volumeCenter[2];


                    double val = getVoxel(pixelCoord);
                
                
                    if (val > maxInt) {
                        maxInt = val;
                    }
                    //0.95 so it can display near to maximum intensity as well.
                    if (val / max > 0.95) {
                        break;
                    }
                }
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxInt/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxInt > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                //smoothens the image when resolution is lowered during Interaction mode
                for (int l = 0; l <res; l++) {
                    for (int m= 0; m <res; m++) {
                        if (m + i < image.getHeight() && m + i >= 0 && l + j < image.getWidth() && l + j >= 0) {
                            image.setRGB(m + i, l + j,pixelColor);
                        }
                    }
                }
            }
        }
    }
    
    private void composite(double[] viewMatrix) {
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        uVec = new double[3];
        vVec = new double[3];
        viewVec = new double[3];
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.normalize(viewVec);
        
        int res = 1;
        if(interactiveMode) {res = 3;}
        // image is square
        imageCenter = image.getWidth() / 2;

        pixelCoord = new double[3];
        volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        voxelColor= new TFColor(); 
        TFColor voxelColorNew = new TFColor();
        

        pixelCoord = new double[3];

        for (int j = 0; j < image.getHeight(); j=j+res) {
            for (int i = 0; i < image.getWidth(); i=i+res) {
                voxelColorNew.set(0, 0, 0, 1);
                int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()); 
                
                for (int step=0;step < diag - 1;step = step + 1) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0]*(step - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1]*(step - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2]*(step - imageCenter) + volumeCenter[2];
                
                double val = getVoxel(pixelCoord);

                // retreiving colour
                voxelColor = tFunc.getColor((int)val);

                    
                // compositing by color multiplication
                voxelColorNew.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * voxelColorNew.r;
                voxelColorNew.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * voxelColorNew.g;
                voxelColorNew.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * voxelColorNew.b;
                }
 
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColorNew.a <= 1.0 ? (int) Math.floor(voxelColorNew.a * 255) : 255;
                int c_red = voxelColorNew.r <= 1.0 ? (int) Math.floor(voxelColorNew.r * 255) : 255;
                int c_green = voxelColorNew.g <= 1.0 ? (int) Math.floor(voxelColorNew.g * 255) : 255;
                int c_blue = voxelColorNew.b <= 1.0 ? (int) Math.floor(voxelColorNew.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                //smoothens the image when resolution is lowered during Interaction mode
                for (int l = 0; l <res; l++) {
                    for (int m= 0; m <res; m++) {
                        if (m + i < image.getHeight() && m + i >= 0 && l + j < image.getWidth() && l + j >= 0) {
                            image.setRGB(m + i, l + j,pixelColor);
                        }
                    }
                }
            }
        }
    }
    
    void twoDTrans (double[] viewMatrix){
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        int res = 1;
        if (interactiveMode) res = 3;
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        viewVec = new double[3];
        uVec = new double[3];
        vVec = new double[3];  
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        imageCenter = image.getWidth() / 2;
        pixelCoord = new double[3];
        volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        max = volume.getMaximum();
        int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ());

        for (int j = 0; j < image.getHeight(); j+=res) {
            for (int i = 0; i < image.getWidth(); i+=res) {
                voxelColor = new TFColor(0,0,0,1);
                
                for(int k = 0; k < diag - 1; k++) {
                    // Calculate the coordinates of the pixel in the data
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (k - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (k - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (k - imageCenter) + volumeCenter[2];
                  
                    int val = (int) getVoxel(pixelCoord);
                    
                    if (val == 0) {
                        continue;
                    }
                    float baseIntensity = tfEditor2D.triangleWidget.baseIntensity;
                    TFColor color = tfEditor2D.triangleWidget.color;
                    double radius = tfEditor2D.triangleWidget.radius;
                    double a;
                    
                    VoxelGradient grad = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);
                    
                    //Kniss approach. Uncomment and set the range to see the effects.
                    /*if (gradient.mag < 0 || gradient.mag > 60)
                        alpha =0;
                    else*/ if (val == baseIntensity && grad.mag == 0){
                        a = 1;
                    } else if (grad.mag > 0 && ((val - (radius * grad.mag)) <= baseIntensity && baseIntensity <= (val + (radius * grad.mag)))) {
                            a = color.a * (1 - (1 / radius) * Math.abs(((baseIntensity - val) / grad.mag))); //Levoy
                        } else {
                        a = 0;
                    }
                    
                    if (shading) {
                        color = getShading(pixelCoord, color,viewVec);
                    }
                    voxelColor.r = (color.r * a) + (voxelColor.r * (1 - a));
                    voxelColor.g = (color.g * a) + (voxelColor.g * (1 - a));
                    voxelColor.b = (color.b * a) + (voxelColor.b * (1 - a));   
                }
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                
                //smoothens the image when resolution is lowered during Interaction mode
                for (int l = 0; l <res; l++) {
                    for (int m= 0; m <res; m++) {
                        if (m + i < image.getHeight() && m + i >= 0 && l + j < image.getWidth() && l + j >= 0) {
                            image.setRGB(m + i, l + j,pixelColor);
                        }
                    }
                }
            }
        }
    }
    
    private TFColor getShading(double[] coord, TFColor voxelColor, double[] viewVec) {
        if (coord[0] < 0 || coord[0] >= volume.getDimX()
                || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            voxelColor.set(0, 0, 0, voxelColor.a);
            return voxelColor;
        }
        TFColor diffuse = voxelColor;
        VoxelGradient vg = gradients.getGradient((int) Math.floor(coord[0]), (int) Math.floor(coord[1]), (int) Math.floor(coord[2]));
        TFColor result = tFunc.getC(viewVec, vg.getNormal(), viewVec, diffuse);
        return result;
        
    }
    

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        //slicer(viewMatrix);    
        switch(opt){
            case 1: slicer(viewMatrix); break;
            case 2: mip(viewMatrix); break;
            case 3: composite(viewMatrix); break;
            case 4: twoDTrans(viewMatrix);break;
            default:
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
    
    

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    
    public void setOpt(int opt){
        this.opt=opt;
        this.changed();
    }
    
    public void setShading(){
        if(shading == true) {
            shading=false;
            this.changed();
        }
        else {
            shading = true;
            this.changed();
        }
    }
}
