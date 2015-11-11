/* ===========================================================
 * JFreeChart : a free chart library for the Java(tm) platform
 * ===========================================================
 *
 * (C) Copyright 2000-2007, by Object Refinery Limited and Contributors.
 *
 * Project Info:  http://www.jfree.org/jfreechart/index.html
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
 * [Java is a trademark or registered trademark of Sun Microsystems, Inc.
 * in the United States and other countries.]
 *
 * -------------------
 * ColorPaintScale.java
 * -------------------
 * (C) Copyright 2006, 2007, by Object Refinery Limited.
 *
 * Original Author:  David Gilbert (for Object Refinery Limited);
 * Contributor(s):   Todd Riley;
 *
 * $Id: ColorPaintScale.java,v 1.1.2.1 2007/01/31 14:15:16 mungady Exp $
 *
 * Changes
 * -------
 * 05-Jul-2006 : Version 1 (DG);
 * 31-Jan-2007 : Renamed min and max to lowerBound and upperBound (DG);
 *
 */

//package org.jfree.chart.renderer;
import org.jfree.chart.renderer.*;

import java.awt.Color;
import java.awt.Paint;
import java.io.Serializable;
import java.util.HashMap;

import org.jfree.util.PublicCloneable;

/**
 * A paint scale that returns colors:
 *  - green to black to red
 *  - black to color
 *  - white to color
 *  - color1 to color2
 *  - rainbow color spectrum
 *
 * @since 1.0.4
 */
public class ColorPaintScale
        implements PaintScale, PublicCloneable, Serializable {

    /** The lower bound. */
    private double lowerBound;

    /** The upper bound. */
    private double upperBound;

    private Color lowerBoundColor;
    private Color upperBoundColor;

    private boolean whiteBackground;
    private boolean rainbow = false;

    HashMap colorHashMap = new HashMap(); // save all new Color objects for re-use


    /**
     * Creates a new <code>ColorPaintScale</code> instance with default values.
     */
    public ColorPaintScale() {
        this(0.0, 1.0, null, null, false, false);
    }

    // default is GreenToBlackToRed
    public ColorPaintScale(double lowerBound, double upperBound) {
        this(lowerBound, upperBound, null, null, false, false);
    }

    public ColorPaintScale(double lowerBound, double upperBound, boolean rainbow) {
        this(lowerBound, upperBound, null, null, false, rainbow);
    }

    // BlackToColor
    public ColorPaintScale(double lowerBound, double upperBound, Color upperBoundColor) {
        this(lowerBound, upperBound, null, upperBoundColor, false, false);
    }

    // BlackToColor or WhiteToColor
    public ColorPaintScale(double lowerBound, double upperBound, Color upperBoundColor, boolean whiteBackground) {
        this(lowerBound, upperBound, null, upperBoundColor, whiteBackground, false);
    }

    // Color1ToColor2
    public ColorPaintScale(double lowerBound, double upperBound, Color lowerBoundColor, Color upperBoundColor) {
        this(lowerBound, upperBound, lowerBoundColor, upperBoundColor, false, false);
    }

    /**
     * Creates a new paint scale for values in the specified range.
     *
     * @param lowerBound  the lower bound.
     * @param upperBound  the upper bound.
     */
    // color scale is BlackToUpperBoundColor
    // black or white background
    public ColorPaintScale(double lowerBound, double upperBound, Color lowerBoundColor, Color upperBoundColor, boolean whiteBackground, boolean rainbow) {
        if (lowerBound >= upperBound) {
            throw new IllegalArgumentException(
                    "Requires lowerBound < upperBound.");
        }
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.upperBoundColor = upperBoundColor;
        this.lowerBoundColor = lowerBoundColor;
        this.whiteBackground = whiteBackground;
        this.rainbow = rainbow;
    }

    /**
     * Returns the lower bound.
     *
     * @return The lower bound.
     */
    public double getLowerBound() {
        return this.lowerBound;
    }

    /**
     * Returns the upper bound.
     *
     * @return The upper bound.
     */
    public double getUpperBound() {
        return this.upperBound;
    }

    /**
     * Returns a paint for the specified value.
     *
     * @param value  the value.
     *
     * @return A paint for the specified value.
     */
    public Paint getPaint(double value) {
        double v = Math.max(value, this.lowerBound);
        v = Math.min(v, this.upperBound);
        //int g = (int) ((value - this.lowerBound) / (this.upperBound - this.lowerBound) * 255.0);
        double g = ((value - this.lowerBound) / (this.upperBound - this.lowerBound));
        if (g > 1.0) {
            g = 1.0;
        }
        if (g < 0.0) {
            g = 0.0;
        }

        int paintRed;
        int paintGreen;
        int paintBlue ;

        if ((this.lowerBoundColor == null) && (this.upperBoundColor != null)) {
            if (whiteBackground) {
                // Paint goes from min=white to max=upperBoundColor
                paintRed   = (int) (((1-g)*Color.white.getRed())   + (g*this.upperBoundColor.getRed()));
                paintGreen = (int) (((1-g)*Color.white.getGreen()) + (g*this.upperBoundColor.getGreen()));
                paintBlue  = (int) (((1-g)*Color.white.getBlue())  + (g*this.upperBoundColor.getBlue())) ;
            }
            else {
                // Paint goes from min=black to max=upperBoundColor
                paintRed   = (int) (g*this.upperBoundColor.getRed());
                paintGreen = (int) (g*this.upperBoundColor.getGreen());
                paintBlue  = (int) (g*this.upperBoundColor.getBlue()) ;
            }
        }
        else if ((this.lowerBoundColor != null) && (this.upperBoundColor != null)) {
            // Paint goes from min=lowerBoundColor to max=upperBoundColor
            paintRed   = (int) (((1-g)*this.lowerBoundColor.getRed())   + (g*this.upperBoundColor.getRed()));
            paintGreen = (int) (((1-g)*this.lowerBoundColor.getGreen()) + (g*this.upperBoundColor.getGreen()));
            paintBlue  = (int) (((1-g)*this.lowerBoundColor.getBlue())  + (g*this.upperBoundColor.getBlue()));

            if (paintRed > 255)   {paintRed = 255;}
            if (paintGreen > 255) {paintGreen = 255;}
            if (paintBlue > 255)  {paintBlue = 255;}
        }
        else if (rainbow == false) { // Paint goes from min=Green to max=Red, avg=black

            //             paintRed = (int)(g*255.0);
            //             paintGreen   = 255 - (int)(g*255.0);
            //             paintBlue  = 0;

            // Paint goes from min=Green to max=Red, avg=black
            if (g > 0.5) {
                paintGreen = 0;
                paintRed = (int) ((g - 0.5)*2.0*255.0);
            }
            else {
                paintRed=0;
                paintGreen = 255 - (int)(g*2.0*255.0);
            }
            paintBlue  = 0;

        }
        else { // rainbow!!
            // there are 5 successive gradients
            if (g <= .2) {
                paintRed = 255 - (int)((g-(0))*5.0*255.0); // ramp down red
                paintGreen = 0;
                paintBlue = 255; // start all blue
            }
            else if (g <= .4) {
                paintRed = 0;
                paintGreen = (int)((g-(.2))*5.0*255.0); // ramp up green
                paintBlue = 255;
            }
            else if (g <= .6) {
                paintRed = 0;
                paintGreen = 255;
                paintBlue = 255 - (int)((g-(.4))*5.0*255.0); // ramp down blue
            }
            else if (g <= .8) {
                paintRed = (int)((g-(.6))*5.0*255.0); // ramp up red;
                paintGreen = 255;
                paintBlue = 0;
            }
            else {
                paintRed = 255;
                paintGreen = 255 - (int)((g-(.8))*5.0*255.0); // ramp down green;
                paintBlue = 0;
            }

            if (paintRed < 0)   {paintRed = 0;}
            if (paintGreen < 0) {paintGreen = 0;}
            if (paintBlue < 0)  {paintBlue = 0;}

            if (paintRed > 255)   {paintRed = 255;}
            if (paintGreen > 255) {paintGreen = 255;}
            if (paintBlue > 255)  {paintBlue = 255;}
        }

        //System.out.println("min="+this.lowerBound+" max="+this.upperBound+" g="+g+" paintRed="+paintRed+" paintGreen="+paintGreen+" paintBlue="+paintBlue+"\n");

        Color aColor = (Color) colorHashMap.get((paintRed*256*256) + (paintGreen*256) + paintBlue);

        if (aColor == null) {
            aColor = new Color(paintRed, paintGreen, paintBlue);
            colorHashMap.put((paintRed*256*256) + (paintGreen*256) + paintBlue, aColor);
        }
        return(aColor);
    }

    /**
     * Tests this <code>ColorPaintScale</code> instance for equality with an
     * arbitrary object.  This method returns <code>true</code> if and only
     * if:
     * <ul>
     * <li><code>obj</code> is not <code>null</code>;</li>
     * <li><code>obj</code> is an instance of <code>ColorPaintScale</code>;</li>
     * </ul>
     *
     * @param obj  the object (<code>null</code> permitted).
     *
     * @return A boolean.
     */
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof ColorPaintScale)) {
            return false;
        }
        ColorPaintScale that = (ColorPaintScale) obj;
        if (this.lowerBound != that.lowerBound) {
            return false;
        }
        if (this.upperBound != that.upperBound) {
            return false;
        }
        return true;
    }

    /**
     * Returns a clone of this <code>ColorPaintScale</code> instance.
     *
     * @return A clone.
     *
     * @throws CloneNotSupportedException if there is a problem cloning this
     *     instance.
     */
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

}
