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
import gui.TransferFunction2DEditor.TriangleWidget;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    public static int renderingMode;

    // ambient reflection coefficient, assuming light source is white
    TFColor SHADING_AMBIENT_COEFF = new TFColor(0.1, 0.1, 0.1, 1.0);
    // diffuse reflection coefficient
    double SHADING_DIFF_COEFF = 0.7;
    // specular reflection coefficient
    double SHADING_SPEC_COEFF = 0.2;
    // exponent used to approximate highligh
    double SHADING_ALPHA = 10;
    
    public boolean shading = false;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
        renderingMode = 0;
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
        tFunc.setTestFunc();

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

    // get voxel intensity
    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] + 1 > volume.getDimX() || coord[1] < 0 || coord[1] + 1 > volume.getDimY()
                || coord[2] < 0 || coord[2] + 1 > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        // speed up
        if(interactiveMode) {
            return volume.getVoxel(x,y,z);
        }
        
        /* trilinear interpolation */
        // calculate alpha, beta and gamma
        double alpha = coord[0] - x;
        double beta = coord[1] - y;
        double gamma = coord[2] - z;

        // intensity values of 8 adjacent pre-defined vetices surrounding the interpolation point
        int sv0 = volume.getVoxel(x, y, z);
        int sv1 = volume.getVoxel(x + 1, y, z);
        int sv2 = volume.getVoxel(x, y + 1, z);
        int sv3 = volume.getVoxel(x + 1, y + 1, z);
        int sv4 = volume.getVoxel(x, y, z + 1);
        int sv5 = volume.getVoxel(x + 1, y, z + 1);
        int sv6 = volume.getVoxel(x, y + 1, z + 1);
        int sv7 = volume.getVoxel(x + 1, y + 1, z + 1);

        // calculate intensity of the given voxel v
        double sv = (1 - alpha) * (1 - beta) * (1 - gamma) * sv0
                + alpha * (1 - beta) * (1 - gamma) * sv1
                + (1 - alpha) * beta * (1 - gamma) * sv2
                + alpha * beta * (1 - gamma) * sv3
                + (1 - alpha) * (1 - beta) * gamma * sv4
                + alpha * (1 - beta) * gamma * sv5
                + (1 - alpha) * beta * gamma * sv6
                + alpha * beta * gamma * sv7;

        return (short) sv;
    }

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void mip(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        // speed up during interactiveMode by decreasing number of computed pixels 
        int resolution;
        if (interactiveMode) {
            resolution = 5;
        } else {
            resolution = 1;
        }
        
        for (int j = 0; j < image.getHeight() - resolution + 1; j += resolution) {
            for (int i = 0; i < image.getWidth() - resolution + 1; i += resolution) {
                
                int maxVal = 0;                        
                
                double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());                       
                double minRange = maxRange * -1;

                // Retrieve the max value along the ray                
                // Interactive mode sample, set sample step = 3, non-interactive mode, step = 1 
                int step = interactiveMode ? 1 : 3;
               
                for (double n = minRange; n < maxRange; n += step) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                                    + viewVec[0] * (n - (maxRange / 2)) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                                    + viewVec[1] * (n - (maxRange / 2)) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                                    + viewVec[2] * (n - (maxRange / 2)) + volumeCenter[2];

                    int val = getVoxel(pixelCoord);

                    if (val > maxVal) {
                        maxVal = val;
                    }
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxVal/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(maxVal);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                        
                for (int ri = 0; ri < resolution; ri++) {
                    for (int rj = 0; rj < resolution; rj++) {
                        if ((i + ri < image.getHeight()) && (j + rj < image.getWidth())) {
                            image.setRGB(ri + i, rj + j, pixelColor);
                        }
                    }
                }                                                         
            }
        }
    }
    
    void compositing(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        // speed up during interactiveMode by decreasing number of computed pixels 
        int resolution;
        if (interactiveMode) {
            resolution = 3;
        } else {
            resolution = 1;
        }
        
        for (int i = 0; i < image.getWidth(); i += resolution) {
            for (int j = 0; j < image.getHeight(); j += resolution) {
               
                TFColor compColor = new TFColor(0,0,0,0);
                double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());
                double minRange = maxRange * -1;

                
                for (double n = minRange; n < maxRange; n += resolution) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * (n - (maxRange / 2)) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * (n - (maxRange / 2)) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * (n - (maxRange / 2)) + volumeCenter[2];

                    int val = getVoxel(pixelCoord);

                    // Apply the transfer function to obtain a color
                    voxelColor = tFunc.getColor(val);

                    compColor.a = voxelColor.a * voxelColor.a + (1 - voxelColor.a) * compColor.a;
                    compColor.r = voxelColor.r * voxelColor.a + (1 - voxelColor.a) * compColor.r;
                    compColor.g = voxelColor.g * voxelColor.a + (1 - voxelColor.a) * compColor.g;
                    compColor.b = voxelColor.b * voxelColor.a + (1 - voxelColor.a) * compColor.b;
                }

                // BufferedImage expects a pixel color packed as ARGB in an int;
                int c_alpha = compColor.a <= 1.0 ? (int) Math.floor(compColor.a * 255) : 255;
                int c_red = compColor.r <= 1.0 ? (int) Math.floor(compColor.r * 255) : 255;
                int c_green = compColor.g <= 1.0 ? (int) Math.floor(compColor.g * 255) : 255;
                int c_blue = compColor.b <= 1.0 ? (int) Math.floor(compColor.b * 255) : 255;

                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                for (int ri = 0; ri < resolution; ri++) {
                    for (int rj = 0; rj < resolution; rj++) {
                        if ((i + ri < image.getHeight()) && (j + rj < image.getWidth())) {
                            image.setRGB(ri + i, rj + j, pixelColor);
                        }
                    }
                }
            }
        }
    }
    
    void transfer(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // parameters to be used to shade
        double [] NormVec = new double[3];
        double dotProductNL = 0;
        double dotProductNH = 0;
        double colorIncrement = 0;
        
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int resolution;
        if (interactiveMode) {
            resolution = 30;
        } else {
            resolution = 1;
        }

        // sample on a plane through all slices of the volume data 
        for (int j = 0; j < image.getHeight(); j += resolution) {
            for (int i = 0; i < image.getWidth(); i += resolution) {
                // Initial color along ray
                voxelColor = new TFColor(0, 0, 0, 1);
                
                double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());
                double minRange = maxRange*-1;
                // find the composite value along the viewing ray
                for (double n = minRange; n < maxRange; n += resolution) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + viewVec[0] * n;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + viewVec[1] * n;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + viewVec[2] * n;

                    int val = getVoxel(pixelCoord);

                    if (val == 0) {
                        continue;
                    }

                    TFColor color = tFunc.getColor(val);

                    TransferFunction2DEditor.TriangleWidget tw = getTF2DPanel().triangleWidget;
                    // user defined parameters
                    short baseIntensity = tw.baseIntensity;
                    double radius = tw.radius;
                    TFColor setColor = tw.color;

                    double minMagnitude = tw.minMagnitude;
                    double maxMagnitude = tw.maxMagnitude;
                    
                    VoxelGradient gradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);
                    
                    color.r = setColor.r;
                    color.g = setColor.g;
                    color.b = setColor.b;
                    
                    // gradient-based opacity weighting
                    if (minMagnitude <= gradient.mag && maxMagnitude >= gradient.mag) {
                         if (val == baseIntensity && gradient.mag == 0){
                             color.a = setColor.a;
                         } else if (gradient.mag > 0 &&
                                 val - (radius * gradient.mag) <= baseIntensity && 
                                 baseIntensity <= val + (radius * gradient.mag)) {
                             color.a = setColor.a * (1 - (1 / radius) * Math.abs((baseIntensity - val) / gradient.mag));
                         } else {
                             color.a = 0;
                         }
                    } else {
                        color.a = 0;
                    }
                    
                    // shading
                    if (shading) {
                        // surface normal at voxel
                        NormVec[0] = gradient.x;
                        NormVec[1] = gradient.y;
                        NormVec[2] = gradient.z;
                        dotProductNL = Math.max(VectorMath.dotproduct(viewVec, NormVec) / (VectorMath.length(NormVec) + 1e-6), 0);
                        dotProductNH = dotProductNL;
                        colorIncrement = SHADING_SPEC_COEFF * Math.pow(dotProductNH, SHADING_ALPHA);
                        color.r = SHADING_AMBIENT_COEFF.r + color.r * (SHADING_DIFF_COEFF * dotProductNL) + colorIncrement;
                        color.g = SHADING_AMBIENT_COEFF.g + color.g * (SHADING_DIFF_COEFF * dotProductNL) + colorIncrement;
                        color.b = SHADING_AMBIENT_COEFF.b + color.b * (SHADING_DIFF_COEFF * dotProductNL) + colorIncrement;
                    }

                    voxelColor.r = (color.r * color.a) + (voxelColor.r * (1 - color.a));
                    voxelColor.g = (color.g * color.a) + (voxelColor.g * (1 - color.a));
                    voxelColor.b = (color.b * color.a) + (voxelColor.b * (1 - color.a));
                    voxelColor.a = (1 - color.a) * voxelColor.a + color.a;
                }

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                
                // Set multiple pixels at lower resolution
                for (int ri = 0; ri < resolution; ri++) {
                    for (int rj = 0; rj < resolution; rj++) {
                        if ((i + ri < image.getHeight()) && (j + rj < image.getWidth())) {
                            image.setRGB(ri + i, rj + j, pixelColor);
                        }
                    }
                }
            }
        }
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

        switch (renderingMode) {
        case 0:
            slicer(viewMatrix);
            break;
        case 1:
            mip(viewMatrix);
            break;
        case 2:
            compositing(viewMatrix);
            break;
        case 3:
            transfer(viewMatrix);
            break;
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
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
