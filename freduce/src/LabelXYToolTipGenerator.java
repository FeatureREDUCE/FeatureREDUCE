/* ===========================================================
 * JFreeChart : a free chart library for the Java(tm) platform
 * ===========================================================
 *
 * (C) Copyright 2000-2008, by Object Refinery Limited and Contributors.
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
 * -------------------------------
 * LabelXYItemLabelGenerator.java
 * -------------------------------
 * (C) Copyright 2002-2008, by Richard Atkinson and Contributors.
 *
 * Original Author:  Richard Atkinson, Todd Riley;
 * Contributor(s):   David Gilbert (for Object Refinery Limited);
 *
 * Changes:
 * --------
 * 05-Aug-2002 : Version 1, contributed by Richard Atkinson (RA);
 * 26-Sep-2002 : Fixed errors reported by Checkstyle (DG);
 * 21-Mar-2003 : Implemented Serializable (DG);
 * 13-Aug-2003 : Implemented Cloneable (DG);
 * 17-Nov-2003 : Implemented PublicCloneable (DG);
 * 25-Feb-2004 : Renamed XYToolTipGenerator --> XYItemLabelGenerator (DG);
 * 02-Feb-2007 : Removed author tags all over JFreeChart sources (DG);
 *
 */

// package org.jfree.chart.labels;

import java.io.Serializable;
import java.util.List;
import java.text.NumberFormat;

import org.jfree.data.xy.XYDataset;
import org.jfree.util.PublicCloneable;
import org.jfree.chart.labels.XYToolTipGenerator;

/**
 * A tool tip generator that stores custom tooltips. The dataset passed into
 * the generateToolTip method is ignored.
 */
public class LabelXYToolTipGenerator implements XYToolTipGenerator,
        Cloneable, PublicCloneable, Serializable {

    /** For serialization. */
    private static final long serialVersionUID = 8636030004670141362L;

    /** Storage for the tooltip lists. */
    private List toolTipSeries = new java.util.ArrayList();

    /**
     * Default constructor.
     */
    public LabelXYToolTipGenerator() {
        super();
    }

    /**
     * Returns the number of tool tip lists stored by the renderer.
     *
     * @return The list count.
     */
    public int getListCount() {
        return this.toolTipSeries.size();
    }

    /**
     * Returns the number of tool tips in a given list.
     *
     * @param series  the series index (zero based).
     *
     * @return The tooltip count.
     */
    public int getToolTipCount(int series) {

        int result = 0;
        List tooltips = (List) this.toolTipSeries.get(series);
        if (tooltips != null) {
            result = tooltips.size();
        }
        return result;
    }

    /**
     * Returns the tool tip text for an item.
     *
     * @param series  the series index.
     * @param item  the item index.
     *
     * @return The tool tip text.
     */
    public String getToolTipText(int series, int item) {

        String result = null;

        if (series < getListCount()) {
            List tooltips = (List) this.toolTipSeries.get(series);
            if (tooltips != null) {
                if (item < tooltips.size()) {
                    // result = (String) tooltips.get(item);
                    result = (String) tooltips.get(tooltips.size()-1 - item);
                }
            }
        }

        return result;
    }

    /**
     * Adds a list of tooltips for a series.
     *
     * @param toolTips  the list of tool tips.
     */
    public void addToolTipSeries(List toolTips) {
        this.toolTipSeries.add(toolTips);
    }

    /**
     * Removes a list of tooltips for a series.
     *
     * @param series  the series index (zero based).
     */
    public void removeToolTipSeries(int series) {
        this.toolTipSeries.remove(series);
    }

    /**
     * Removes all the tooltips series.
     *
     */
    public void clearToolTipSeries() {
        this.toolTipSeries.clear();
    }

    /**
     * Generates a tool tip text item for a particular item within a series.
     *
     * @param data  the dataset (ignored in this implementation).
     * @param series  the series (zero-based index).
     * @param item  the item (zero-based index).
     *
     * @return The tooltip text.
     */
    public String generateToolTip(XYDataset dataset, int series, int item) {

        String seriesLabel = dataset.getSeriesKey(series).toString();
        Number x = dataset.getX(series, item);
        Number y = dataset.getY(series, item);
        String pointLabel = getToolTipText(series, item);

        String toolTip =
            "<html><center>"
            +"<p style='color:green;'><b>"
            +pointLabel
            +"</b></p>"
            // +"<br>"
            +"<p style='color:purple;'>"
            +seriesLabel
            +"</p>"
            // +"<br>"
            +"("
            +format(x.doubleValue(),2,2)
            +", "
            +format(y.doubleValue(),2,2)
            +")</center></html>";

        return(toolTip);

    }

    private String format(double num, int minDecimalPlaces, int maxDecimalPlaces) {
        NumberFormat format = NumberFormat.getInstance();
        format.setMinimumFractionDigits(minDecimalPlaces);
        format.setMaximumFractionDigits(maxDecimalPlaces);
        return(format.format(num));
    }

    /**
     * Returns an independent copy of the generator.
     *
     * @return A clone.
     *
     * @throws CloneNotSupportedException if cloning is not supported.
     */
    public Object clone() throws CloneNotSupportedException {

        LabelXYToolTipGenerator clone
            = (LabelXYToolTipGenerator) super.clone();
        if (this.toolTipSeries != null) {
            clone.toolTipSeries = new java.util.ArrayList(this.toolTipSeries);
        }
        return clone;

    }

    /**
     * Tests if this object is equal to another.
     *
     * @param obj  the other object.
     *
     * @return A boolean.
     */
    public boolean equals(Object obj) {

        if (obj == this) {
            return true;
        }

        if (obj instanceof LabelXYToolTipGenerator) {
            LabelXYToolTipGenerator generator = (LabelXYToolTipGenerator) obj;
            boolean result = true;
            for (int series = 0; series < getListCount(); series++) {
                for (int item = 0; item < getToolTipCount(series); item++) {
                    String t1 = getToolTipText(series, item);
                    String t2 = generator.getToolTipText(series, item);
                    if (t1 != null) {
                        result = result && t1.equals(t2);
                    }
                    else {
                        result = result && (t2 == null);
                    }
                }
            }
            return result;
        }

        return false;

    }

}
