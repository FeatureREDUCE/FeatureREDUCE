/* ===========================================================
 * BioJavaChart
 * ===========================================================
 *
 *
 * Original Author: Todd Riley
 * Contributor(s):  Matthias Rose (Ablay & Fodi GmbH, Germany);
 *		    Eduardo Ramalho;
 *		    David Gilbert (for Object Refinery Limited);
 *
 * This code uses code from PanScrollZoomDemo.java
 *
 *
 */

import java.awt.Dimension;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.BasicStroke;

import java.awt.event.WindowListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Line2D;

import java.io.IOException;
import java.text.NumberFormat;

import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.event.ChartChangeEvent;
import org.jfree.chart.event.ChartChangeListener;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.ValueAxisPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.GrayPaintScale;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.Range;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYShapeAnnotation;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.xy.XYBarDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.general.SeriesException;

import org.jfree.chart.plot.ContourPlot;
import org.jfree.data.contour.ContourDataset;
import org.jfree.data.contour.DefaultContourDataset;
import org.jfree.data.contour.NonGridContourDataset;
import org.jfree.chart.axis.ColorBar;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.axis.Axis;
import org.jfree.chart.axis.SymbolAxis;
import org.jfree.data.statistics.DefaultStatisticalCategoryDataset;
import org.jfree.data.statistics.StatisticalCategoryDataset;
import org.jfree.data.statistics.HistogramType;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.data.xy.DefaultXYZDataset;
import java.awt.GradientPaint;
import java.awt.Color;

import org.jfree.data.category.CategoryDataset;
import org.jfree.data.general.DatasetUtilities;

import javax.swing.AbstractButton;
import javax.swing.BoundedRangeModel;
import javax.swing.ButtonGroup;
import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.JCheckBox;
import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JScrollBar;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.JFrame;
import javax.swing.JSlider;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import javax.swing.BoxLayout;
import javax.swing.Box;
import java.awt.Component;
import java.io.*;

import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.ItemLabelAnchor;
import org.jfree.chart.labels.ItemLabelPosition;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.StatisticalBarRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.ui.TextAnchor;
import org.jfree.chart.urls.StandardCategoryURLGenerator;
import org.jfree.chart.entity.StandardEntityCollection;
import org.jfree.chart.ChartUtilities;

import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RectangleInsets;

import java.util.LinkedHashMap;
import java.util.List;

import org.biojava.bio.mydp.*;
import org.biojava.bio.mydp.FileTools;
import org.biojava.bio.mydp.StringTools;



/**
 * A demo of the {@link HistogramDataset} class.
 */
public class BioJavaChart extends JPanel
    implements ActionListener,
               ChangeListener,
               ChartChangeListener,
               MouseListener,
               MouseMotionListener {

    /** The panel that displays the chart. */
    public ChartPanel chartPanel;

    public JFreeChart chart = null;

    private JSlider binsSlider = null;

    /** The scroll factor. */
    private double scrollFactor = 1000;

    /** The scroll bar. */
    private JScrollBar scrollBar;

    /** The starting point for panning. */
    private Point2D panStartPoint;

    /** The min/max values for the primary axis. */
    private double[] primYMinMax = new double[2];

    /** The min/max values for the secondary axis. */
    private double[] secondYMinMax = new double[2];

    /** Action command for the 'Pan' button. */
    private static final String ACTION_CMD_PAN = "pan";

    /** Action command for the zoom box button. */
    private static final String ACTION_CMD_ZOOM_BOX = "zoomBox";

    /** Action command for the zoom fit button. */
    private static final String ACTION_CMD_ZOOM_TO_FIT = "zoomFit";

    /** Action command for the '+' button. */
    private static final String ACTION_CMD_ZOOM_IN = "zoomIn";

    /** Action command for the '-' button. */
    private static final String ACTION_CMD_ZOOM_OUT = "zoomOut";

    private static final String ACTION_CMD_X_LINEAR = "xLinear";
    private static final String ACTION_CMD_X_LOG = "xLog";
    private static final String ACTION_CMD_Y_LINEAR = "yLinear";
    private static final String ACTION_CMD_Y_LOG = "yLog";
    private static final String ACTION_CMD_UNITY_LINE = "unityLine";

    /** The zoom factor. */
    private static final double ZOOM_FACTOR = 0.8;

    /** The toolbar. */
    private JToolBar toolBar;

    /** The zoom button. */
    private AbstractButton zoomButton;

    /** The pan button. */
    private AbstractButton panButton;

    /** The zoom in button. */
    private AbstractButton zoomInButton;

    /** The zoom out button. */
    private AbstractButton zoomOutButton;

    /** The fit button. */
    private AbstractButton fitButton;

    // private JRadioButton linearXButton = null;
    // private JRadioButton logXButton = null;
    // private JRadioButton linearYButton = null;
    // private JRadioButton logYButton = null;
    // private JCheckBox unityLineCheckBox = null;

    private AbstractButton linearXButton = null;
    private AbstractButton logXButton = null;
    private AbstractButton linearYButton = null;
    private AbstractButton logYButton = null;
    private AbstractButton unityLineCheckBox = null;

    int maxBins = 1000;

    boolean includeLogLinear = false;

    // Histogram globals
    double[][] histValues = null;
    String[] histDatasetLabels = null;
    IntervalXYDataset histDataset = null;
    boolean numBinsChanged = false;

    // ROC curve globals
    String[] rocDatasetLabels = null;
    String rocPath = null;
    String[] rocFileNames = null;
    XYSeriesCollection rocDataset = null;

    Output out = new Output(Output.Destination.STDOUT, (byte)25, false);

    public BioJavaChart() {
    }

    /**
     * Creates a new demo.
     *
     * @param title  the frame title.
     */
//     public BioJavaChart(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String datasetLabel, boolean legend, double[] values, int bins) {
//         initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabel, legend, values, bins);
//     }

    public BioJavaChart(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String[] datasetLabels, boolean legend, String path, String[] fileNames, int bins) {
        initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, path, fileNames, bins);
    }

    public BioJavaChart(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String[] datasetLabels, boolean legend, double[][] values, int bins) {
        initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, values, bins);
    }

    public BioJavaChart(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String[] datasetLabels, boolean legend, String[][] values, int bins) {
        initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, values, bins);
    }

//     public void initialize(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String datasetLabel, boolean legend,  double[] values, int bins) {
//         //super(title);
//         //chart = null;

//         if (chartType.equalsIgnoreCase("histogram")) {
//             if (bins == -1) {
//                 XYBarDataset dataset = createXYBarDataset(datasetLabel, values);
//                 chart = createXYBarChart(title, xAxisLabel, yAxisLabel, legend, dataset);
//             }
//             else {
//                 IntervalXYDataset dataset = createHistDataset(datasetLabel, values, bins);
//                 chart = createHistChart(title, xAxisLabel, yAxisLabel, legend, dataset);
//             }
//         }

//         if (chart != null) {
//             createPanScrollZoomPanel(chart);
//         }

//         initialize();

//     }

    public void initialize(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String[] datasetLabels, boolean legend,  String path, String[] fileNames, int bins) {
        //super(title);
        //chart = null;

        boolean includeBinsSlider = false;
        preInitialize();

        if (chartType.equalsIgnoreCase("roc")) {

            includeBinsSlider = true;
            if (title == null) {
                title = "Receiver Operating Characteristic (ROC) Curve";
            }
            if (xAxisLabel == null) {
                xAxisLabel = "False Positive Rate";
            }
            if (yAxisLabel == null) {
                yAxisLabel = "True Positive Rate";
            }

            //XYSeriesCollection rocDataset = new XYSeriesCollection();
            rocDataset = new XYSeriesCollection();

            RocData[] rocDataArray = new RocData[fileNames.length];
            for (int i=0; i < fileNames.length; i++) {
                if (path != null) {
                    rocDataArray[i] = new RocData(datasetLabels[i], path + File.separator + fileNames[i], bins);
                }
                else {
                    rocDataArray[i] = new RocData(datasetLabels[i], fileNames[i], bins);
                }
                String aucString = ": A.U.C. = "+format(rocDataArray[i].auc, 3, 3);

                XYSeries series = null;
                if (datasetLabels != null) {
                    series = createXYSeries(datasetLabels[i]+aucString, rocDataArray[i].rocArray[0], rocDataArray[i].rocArray[1]);
                }
                else {
                    series = createXYSeries(aucString, rocDataArray[i].rocArray[0], rocDataArray[i].rocArray[1]);
                }
                rocDataset.addSeries(series);

                //createPanScrollZoomPanel(chart);
            }

            chart = createLineChart(title, xAxisLabel, yAxisLabel, legend, rocDataset);

            this.rocDatasetLabels = datasetLabels;
            this.rocPath = path;
            this.rocFileNames = fileNames;

        }

        if (chart != null) {
            createPanScrollZoomPanel(chart, includeBinsSlider);
        }

        postInitialize();

    }

    public void initialize(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String[] datasetLabels, boolean legend, double[][] values, int bins) {
        //JFreeChart chart = null;
        //JPanel panel = null;

        boolean includeBinsSlider = false;
        preInitialize();

        if (chartType.equalsIgnoreCase("scatter") || chartType.equalsIgnoreCase("lines")) {
            includeLogLinear = true;
            XYSeriesCollection dataset = createXYSeriesDataset(datasetLabels, values);

            if (chartType.equalsIgnoreCase("scatter")) {
                chart = createScatterChart(title, xAxisLabel, yAxisLabel, legend, dataset);
            }
            else if (chartType.equalsIgnoreCase("lines")) {
                chart = createLineChart(title, xAxisLabel, yAxisLabel, legend, dataset);
            }
        }
        else if (chartType.equalsIgnoreCase("histogram")) {
            if (bins == -1) {
                XYBarDataset dataset = createXYBarDataset(datasetLabels, values, false);
                chart = createXYBarChart(title, xAxisLabel, yAxisLabel, legend, dataset);
            }
            else {
                includeBinsSlider = true;

                histDataset = createHistDataset(datasetLabels, values, bins);
                chart = createHistChart(title, xAxisLabel, yAxisLabel, legend, histDataset);

                // IntervalXYDataset dataset = createHistDataset(datasetLabels, values, bins);
                // chart = createHistChart(title, xAxisLabel, yAxisLabel, legend, dataset);

                this.histValues = values;
                this.histDatasetLabels = datasetLabels;
                //histDataset = dataset;
            }
        }

        if (chart != null) {
            createPanScrollZoomPanel(chart, includeBinsSlider);
        }

        postInitialize();

    }

    public void initialize(String chartType, String title, String xAxisLabel, String[] xAxisStrings, String yAxisLabel, String[] yAxisStrings, String zAxisLabel, String[] zAxisStrings, String[] datasetLabels, boolean legend, double[][] values, int bins) {

        if (chartType.equalsIgnoreCase("heatmap")) {

            preInitialize();

            XYZDataset dataset = null;
            //ContourDataset dataset = null;
            if (datasetLabels != null) {
                dataset = createXYZDataset(datasetLabels[0], values);
                //dataset = createNonGridContourDataset(datasetLabels[0], values[0], values[1], values[2]);
            }
            else {
                dataset = createXYZDataset(null, values);
                //dataset = createNonGridContourDataset(null, values[0], values[1], values[2]);
            }
            //chart = createContourChart(title, xAxisLabel, yAxisLabel, zAxisLabel, legend, dataset);

            double minZ = MathTools.min(values[2]);
            double maxZ = MathTools.max(values[2]);

//            chart = createXYBlockChart(title, xAxisLabel, yAxisLabel, zAxisLabel, legend, dataset, minZ, maxZ);
            chart = createXYBlockChart(title, xAxisLabel, xAxisStrings, yAxisLabel, yAxisStrings, zAxisLabel, zAxisStrings, legend, dataset, minZ, maxZ);
            createPanScrollZoomPanel(chart);

            postInitialize();

        }
    }

    public void initialize(String chartType, String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, String[] datasetLabels, boolean legend, String[][] values, int bins) {
        //JFreeChart chart = null;
        //JPanel panel = null;

        preInitialize();

        if (chartType.equalsIgnoreCase("bar")) {
            CategoryDataset dataset = null;

            // If there are no St. Dev. values
            if (values.length == 2) {
                if (datasetLabels != null) {
                    dataset = createCategoryDataset(datasetLabels[0], values[0], StringTools.toDoubles(values[1]));
                }
                else {
                    dataset = createCategoryDataset(null, values[0], StringTools.toDoubles(values[1]));
                }
            }
            else {
                if (datasetLabels != null) {
                    dataset = createCategoryDataset(datasetLabels[0], values[0], StringTools.toDoubles(values[1]), StringTools.toDoubles(values[2]));
                }
                else {
                    dataset = createCategoryDataset(null, values[0], StringTools.toDoubles(values[1]), StringTools.toDoubles(values[2]));
                }
            }

            chart = createBarChart(title, xAxisLabel, yAxisLabel, legend, dataset);
            createDefaultChartPanel(chart);
        }
        else {
            initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, StringTools.toDoubleMatrix(values), bins);
        }

        postInitialize();
    }

    private void preInitialize() {
        ChartFactory.setChartTheme(StandardChartTheme.createLegacyTheme());
        BarRenderer.setDefaultBarPainter(new StandardBarPainter());
        BarRenderer.setDefaultShadowsVisible(false);
        XYBarRenderer.setDefaultBarPainter(new StandardXYBarPainter());
        XYBarRenderer.setDefaultShadowsVisible(false);
    }

    // This initialize() routine should always be called by the others!!!!!!!!!!!!!!!
    private void postInitialize() {

        //Set background to white so that the PNG graphics look better for publications
        if (chart != null) {
            chart.setBackgroundPaint(Color.white);
        }
        else {
            out.println("In initialize(): chart is null!\n");
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    private CategoryDataset createCategoryDataset(String datasetLabel, String[] columnKeys, double[] data) {
        String[] rowKeys = new String[1];
        if (datasetLabel != null) {
            rowKeys[0] = datasetLabel;
        }
        else {
            rowKeys[0] = "Dataset 1";
        }

        double[][] newData = new double[1][data.length];
        newData[0] = data;
        CategoryDataset dataset = DatasetUtilities.createCategoryDataset(rowKeys, columnKeys, newData);
        return dataset;
    }

    private CategoryDataset createCategoryDataset(String[] datasetLabels, String[] columnKeys, double[][] data) {
        CategoryDataset dataset = DatasetUtilities.createCategoryDataset(datasetLabels, columnKeys, data);
        return dataset;
    }

    private CategoryDataset createCategoryDataset(String datasetLabel, String[] columnKeys, double[] means, double[] stdevs) {
        String[] rowKeys = new String[1];
        if (datasetLabel != null) {
            rowKeys[0] = datasetLabel;
        }
        else {
            rowKeys[0] = "Dataset 1";
        }

        DefaultStatisticalCategoryDataset dataset = new DefaultStatisticalCategoryDataset();
        for (int j=0; j < means.length; j++) {
            dataset.add(means[j], stdevs[j], rowKeys[0], columnKeys[j]);
        }
        return dataset;
    }


    // Modified source from 1.0.2 ChartFactory.java
    // modification uses StatisticalBarRenderer for a StatisticalCategoryDataset
    public static JFreeChart createBarChart(String title,
        String categoryAxisLabel,
        String valueAxisLabel,
        CategoryDataset dataset,
        PlotOrientation orientation,
        boolean legend,
        boolean tooltips,
        boolean urls) {

        if (orientation == null) {
            throw new IllegalArgumentException("Null 'orientation' argument.");
        }
        CategoryAxis categoryAxis = new CategoryAxis(categoryAxisLabel);
        ValueAxis valueAxis = new NumberAxis(valueAxisLabel);

        BarRenderer renderer;
        if (dataset instanceof StatisticalCategoryDataset) {
            renderer = new StatisticalBarRenderer();
        }
        else {
            renderer = new BarRenderer();
        }

        if (orientation == PlotOrientation.HORIZONTAL) {
            ItemLabelPosition position1 = new ItemLabelPosition(
                ItemLabelAnchor.OUTSIDE3, TextAnchor.CENTER_LEFT);
            renderer.setPositiveItemLabelPosition(position1);
            ItemLabelPosition position2 = new ItemLabelPosition(
                ItemLabelAnchor.OUTSIDE9, TextAnchor.CENTER_RIGHT);
            renderer.setNegativeItemLabelPosition(position2);
        }
        else if (orientation == PlotOrientation.VERTICAL) {
            ItemLabelPosition position1 = new ItemLabelPosition(
                ItemLabelAnchor.OUTSIDE12, TextAnchor.BOTTOM_CENTER);
            renderer.setPositiveItemLabelPosition(position1);
            ItemLabelPosition position2 = new ItemLabelPosition(
                ItemLabelAnchor.OUTSIDE6, TextAnchor.TOP_CENTER);
            renderer.setNegativeItemLabelPosition(position2);
        }
        if (tooltips) {
            renderer.setBaseToolTipGenerator(
                new StandardCategoryToolTipGenerator());
        }
        if (urls) {
            renderer.setBaseItemURLGenerator(
                new StandardCategoryURLGenerator());
        }

        CategoryPlot plot = new CategoryPlot(dataset, categoryAxis, valueAxis,
            renderer);
        plot.setOrientation(orientation);
        JFreeChart chart = new JFreeChart(title, JFreeChart.DEFAULT_TITLE_FONT,
            plot, legend);

        return chart;

    }

    private JFreeChart createBarChart(String title, String xAxisLabel, String yAxisLabel, boolean legend, CategoryDataset dataset) {
        //JFreeChart chart = ChartFactory.createBarChart(
        JFreeChart chart = createBarChart(
            title,
            xAxisLabel,
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            legend,
            true,
            false
                                          );
        //chart.getXYPlot().setForegroundAlpha(0.75f);
        return chart;
    }
    /////////////////////////////////////////////////////////////////////////////

    /**
     * Creates a sample {@link HistogramDataset}.
     *
     * @return the dataset.
     */
    //     private XYBarDataset createXYBarDataset() {
    //         HistogramDataset dataset = new HistogramDataset();
    //         double[] values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    //         dataset.addSeries("H1", values, 10, 0.0, 10.0);
    //         return (new XYBarDataset(dataset, 1));
    //     }

    public void setXYDataset(String[] datasetLabels, double[][] values) {

        XYSeriesCollection xYSeriesCollection = createXYSeriesDataset(datasetLabels, values);

        chart.getXYPlot().setDataset(xYSeriesCollection);
    }

    public void addXYDataset(String[] datasetLabels, double[][] values) {
        addXYSeriesDataset(datasetLabels, values);
    }

    public void removeXYDatasets() {
        XYSeriesCollection dataset = (XYSeriesCollection) chart.getXYPlot().getDataset();

        dataset.removeAllSeries();
    }

    public void removeXYDataset(int series) {
        XYSeriesCollection dataset = (XYSeriesCollection) chart.getXYPlot().getDataset();

        dataset.removeSeries(series);
    }

    public void addXYToolTips(List<String> toolTipsStrings) {
        XYToolTipGenerator toolTipsGenerator = (XYToolTipGenerator) chart.getXYPlot().getRenderer().getBaseToolTipGenerator();

        if (((toolTipsGenerator == null)) || !(toolTipsGenerator instanceof LabelXYToolTipGenerator)) {
            toolTipsGenerator = new LabelXYToolTipGenerator();
            chart.getXYPlot().getRenderer().setBaseToolTipGenerator(toolTipsGenerator);
        }

        ((LabelXYToolTipGenerator)toolTipsGenerator).addToolTipSeries(toolTipsStrings);
    }

    public void clearXYToolTips() {
        XYToolTipGenerator toolTipsGenerator = (XYToolTipGenerator) chart.getXYPlot().getRenderer().getBaseToolTipGenerator();

        if ((toolTipsGenerator != null) && ((toolTipsGenerator instanceof LabelXYToolTipGenerator))) {
            ((LabelXYToolTipGenerator)toolTipsGenerator).clearToolTipSeries();
        }
    }

    public void removeXYToolTips(int series) {
        XYToolTipGenerator toolTipsGenerator = (XYToolTipGenerator) chart.getXYPlot().getRenderer().getBaseToolTipGenerator();

        if ((toolTipsGenerator != null) && ((toolTipsGenerator instanceof LabelXYToolTipGenerator))) {
            ((LabelXYToolTipGenerator)toolTipsGenerator).removeToolTipSeries(series);
        }
    }

    private XYSeriesCollection createXYSeriesDataset(String[] datasetLabels, double[][] values) {
        XYSeriesCollection dataset = new XYSeriesCollection();

        // There is 1 less dataset label than value columns
        //
        // # dataset1 dataset2   ( 2 labels )
        // x y1 y2               ( 3 values )
        //
        if (values != null) {
            for (int i = 1; i < values.length; i++) {
                XYSeries series = null;
                if (datasetLabels != null) {
                    series = createXYSeries(datasetLabels[i-1], values[0], values[i]);
                }
                else {
                    series = createXYSeries(null, values[0], values[i]);
                }
                dataset.addSeries(series);
            }

        }
        return (dataset);
    }

    private void addXYSeriesDataset(String[] datasetLabels, double[][] values) {
        XYSeriesCollection dataset = (XYSeriesCollection) chart.getXYPlot().getDataset();

        // There is 1 less dataset label than value columns
        //
        // # dataset1 dataset2   ( 2 labels )
        // x y1 y2               ( 3 values )
        //
        if (values != null) {
            for (int i = 1; i < values.length; i++) {
                XYSeries series = null;
                if (datasetLabels != null) {
                    series = createXYSeries(datasetLabels[i-1], values[0], values[i]);
                }
                else {
                    series = createXYSeries(null, values[0], values[i]);
                }
                dataset.addSeries(series);
            }

        }
        return;
    }

    private XYSeries createXYSeries(String datasetLabel, double[] xValuesArray, double[] yValuesArray) {
        XYSeries series = null;

        if (datasetLabel == null) {
            datasetLabel = new String("Dataset 1");
        }
        series = new XYSeries(datasetLabel);

        for (int i = 0; i < xValuesArray.length; i++) {
            try {
                series.add(xValuesArray[i], yValuesArray[i]);
                //System.err.println("\tAdding Series "+i);
            }
            catch (SeriesException e) {
                System.err.println("Error adding to series");
            }
        }
        return(series);
    }

    private JFreeChart createLineChart(String title, String xAxisLabel, String yAxisLabel, boolean legend, XYSeriesCollection dataset) {
        JFreeChart chart = ChartFactory.createXYLineChart(
            title,
            xAxisLabel,
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            legend,
            true,
            false
                                                          );

        XYPlot plot = chart.getXYPlot();
//         plot.setBackgroundPaint(Color.lightGray);
//         plot.setDomainGridlinePaint(Color.white);
//         plot.setRangeGridlinePaint(Color.white);

        // set the domain axis to display integers only...
        final NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        double length = domainAxis.getRange().getLength();
        if (length > 4) {
            domainAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        }

        // set the range axis to display integers only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        length = rangeAxis.getRange().getLength();
        if (length > 4) {
            rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        }

        //chart.getXYPlot().setForegroundAlpha(0.75f);
        return chart;
    }

    private JFreeChart createScatterChart(String title, String xAxisLabel, String yAxisLabel, boolean legend, XYSeriesCollection dataset) {
        JFreeChart chart = ChartFactory.createScatterPlot(
            title,
            xAxisLabel,
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            legend,
            true,
            false
                                                          );
        XYPlot plot = chart.getXYPlot();
//         plot.setBackgroundPaint(Color.lightGray);
//         plot.setDomainGridlinePaint(Color.white);
//         plot.setRangeGridlinePaint(Color.white);

        // set the domain axis to display integers only...
        final NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        double length = domainAxis.getRange().getLength();
        if (length > 4) {
            domainAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        }

        // set the range axis to display integers only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        length = rangeAxis.getRange().getLength();
        if (length > 4) {
            rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        }

        //chart.getXYPlot().setForegroundAlpha(0.75f);
        return chart;
    }

    /////////////////////////////////////////////////////////////////////////////

    private IntervalXYDataset createHistDataset(String[] datasetLabels, double[][] values, int bins) {
        HistogramDataset dataset = new HistogramDataset();
        if (datasetLabels == null) {
            datasetLabels = new String[values.length];
        }

        for (int i=0; i < values.length; i++) {
            if (datasetLabels[i] == null) {
                datasetLabels[i] = new String("Dataset "+(i+1));
            }
            dataset.addSeries(datasetLabels[i], values[i], bins);
        }
        dataset.setType(HistogramType.RELATIVE_FREQUENCY);
        return dataset;
    }

    private JFreeChart createHistChart(String title, String xAxisLabel, String yAxisLabel, boolean legend, IntervalXYDataset dataset) {
        JFreeChart chart = ChartFactory.createHistogram(
            title,
            xAxisLabel,
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            legend,
            true,
            false
                                                        );
        //chart.getXYPlot().setForegroundAlpha(0.75f);
        return chart;
    }
    /////////////////////////////////////////////////////////////////////////////



    // ****************************************************************************
    // * JFREECHART DEVELOPER GUIDE                                               *
    // * The JFreeChart Developer Guide, written by David Gilbert, is available   *
    // * to purchase from Object Refinery Limited:                                *
    // *                                                                          *
    // * http://www.object-refinery.com/jfreechart/guide.html                     *
    // *                                                                          *
    // * Sales are used to provide funding for the JFreeChart project - please    *
    // * support us so that we can continue developing free software.             *
    // ****************************************************************************


    private XYBarDataset createXYBarDataset(String datasetLabels[], double[][] values, boolean freqs) {

        if (datasetLabels == null) {
            datasetLabels = new String[values.length];
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        XYBarDataset barDataset = new XYBarDataset(dataset,1);
        int index;

        //out.println("Length = "+values.length);

        for (int i=0; i < values.length; i++) {

            if (datasetLabels[i] == null) {
                datasetLabels[i] = new String("Dataset "+(i+1));
            }

            XYSeries series = new XYSeries(datasetLabels[i]);

            for (int j=0; j < values[i].length; j++) {
                try {
                    //out.println("Adding Value = "+values[i]);
                    index = series.indexOf(values[i][j]);
                    if (index > -1) {
                        series.updateByIndex(index, series.getY(index).doubleValue() + 1.0 );
                    }
                    else {
                        series.add(values[i][j], 1.0);
                    }
                    //System.err.println("\tAdding Series "+i);
                }
                catch (SeriesException e) {
                    System.err.println("Error adding to series");
                }
            }

            //Go through the counts and make them frequencies
            if (freqs) {
                int itemCount = series.getItemCount();
                for (int j = 0; j < itemCount; j++) {
                    series.updateByIndex(j, series.getY(j).doubleValue() / values[i].length);
                }
            }

            dataset.addSeries(series);

        }

        //out.println("Series Item count = "+series.getItemCount());
        //out.println("AllowDups = "+series.getAllowDuplicateXValues());

        //return (new XYBarDataset(dataset,1));
        return (barDataset);
    }


    /**
     * Creates a chart.
     *
     * @param dataset  a dataset.
     *
     * @return The chart.
     */
    private JFreeChart createXYBarChart(String title, String xAxisLabel, String yAxisLabel, boolean legend, XYBarDataset dataset) {
        JFreeChart chart = ChartFactory.createXYBarChart(
            title,
            xAxisLabel,
            false, //xAxis dates?
            yAxisLabel,
            dataset,
            PlotOrientation.VERTICAL,
            legend,
            true, // tooltips
            false);

        chart.getXYPlot().setForegroundAlpha(0.75f);
        //((XYBarRenderer)chart.getXYPlot().getRenderer()).setDrawBarOutline(true);
        //((XYBarRenderer)chart.getXYPlot().getRenderer()).setDefaultShadowsVisible(false);

        return chart;
    }


    private XYZDataset createXYZDataset(String datasetLabel, double[][] xyzValuesMatrix) {

        DefaultXYZDataset dataset = new DefaultXYZDataset();
        dataset.addSeries(datasetLabel, xyzValuesMatrix);
        return(dataset);

    }


    private XYZDataset createXYZDataset(String datasetLabel, double[] xValuesArray, double[] yValuesArray, double[] zValuesArray) {

        double[][] data = new double[][] {xValuesArray, yValuesArray, zValuesArray};
        return(createXYZDataset(datasetLabel, data));
        // DefaultXYZDataset dataset = new DefaultXYZDataset();
        // dataset.addSeries(datasetLabel, data);
        // return(dataset);

    }

    private ContourDataset createNonGridContourDataset(String datasetLabel, double[] xValuesArray, double[] yValuesArray, double[] zValuesArray) {

        Double[] xDoubles = new Double[xValuesArray.length];
        Double[] yDoubles = new Double[yValuesArray.length];
        Double[] zDoubles = new Double[zValuesArray.length];

        for (int i = 0; i < xValuesArray.length; i++) {
            xDoubles[i] = new Double(xValuesArray[i]);
            yDoubles[i] = new Double(yValuesArray[i]);
            zDoubles[i] = new Double(zValuesArray[i]);
            //zDoubles[i] = new Double(xValuesArray[i]+yValuesArray[i]*1000);
        }

        NonGridContourDataset dataset = new NonGridContourDataset(datasetLabel, xDoubles, yDoubles, zDoubles);

        //NonGridContourDataset dataset = new NonGridContourDataset("Contour Plot", xDoubles, yDoubles, zDoubles, 100, 100, 100);
        //DefaultContourDataset dataset = new DefaultContourDataset("Contour Plot", xDoubles, yDoubles, zDoubles);

        return (dataset);
    }


    /**
     * Creates a chart.
     *
     * @param dataset  the dataset.
     *
     * @return A chart.
     */
    private JFreeChart createXYBlockChart(String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, boolean legend, XYZDataset dataset, double minZ, double maxZ) {
        return(createXYBlockChart(title, xAxisLabel, null, yAxisLabel, null, zAxisLabel, null, legend, dataset, minZ, maxZ));
    }

    private JFreeChart createXYBlockChart(
        String title,
        String xAxisLabel,
        String[] xAxisStrings,
        String yAxisLabel,
        String[] yAxisStrings,
        String zAxisLabel,
        String[] zAxisStrings,
        boolean legend,
        XYZDataset dataset,
        double minZ,
        double maxZ) {

//         ValueAxis xAxis = null;    /** The x-axis. */
//         NumberAxis yAxis = null;   /** The y-axis. */
//         NumberAxis zAxis = null;   /** The z-axis. */

        ValueAxis xAxis = null;    /** The x-axis. */
        ValueAxis yAxis = null;   /** The y-axis. */
        ValueAxis zAxis = null;   /** The z-axis. */

        ///////////////////////////////////////
        // xAxis
        ///////////////////////////////////////
        if (xAxisStrings == null) {
            xAxis = new NumberAxis(xAxisLabel);
            ((NumberAxis) xAxis).setLowerMargin(0.0);
            ((NumberAxis) xAxis).setUpperMargin(0.0);
            ((NumberAxis) xAxis).setAutoRangeIncludesZero(false);
            ((NumberAxis) xAxis).setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        }
        else {
            xAxis = new SymbolAxis(xAxisLabel, xAxisStrings);
            xAxis.setVerticalTickLabels(true);
        }

        ///////////////////////////////////////
        // yAxis
        ///////////////////////////////////////
        if (yAxisStrings == null) {
            yAxis = new NumberAxis(yAxisLabel);
            ((NumberAxis) yAxis).setLowerMargin(0.0);
            ((NumberAxis) yAxis).setUpperMargin(0.0);
            ((NumberAxis) yAxis).setAutoRangeIncludesZero(false);
            ((NumberAxis) yAxis).setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        }
        else {
            yAxis = new SymbolAxis(yAxisLabel, yAxisStrings);
        }

        ///////////////////////////////////////
        // zAxis
        ///////////////////////////////////////
        if (zAxisStrings == null) {
            zAxis = new NumberAxis(zAxisLabel);

            // if minZ or maxZ are very close to an integer, then set them to that integer
            long roundedMinZ = Math.round(minZ);
            long roundedMaxZ = Math.round(maxZ);
            if (Math.abs(roundedMinZ - minZ) < .001) {
                minZ = roundedMinZ;
            }
            if (Math.abs(roundedMaxZ - maxZ) < .001) {
                maxZ = roundedMaxZ;
            }

        }
        else {
            zAxis = new SymbolAxis(zAxisLabel, zAxisStrings);
        }
        ((NumberAxis) zAxis).setRange(minZ, maxZ);

        XYBlockRenderer renderer = new XYBlockRenderer();

        //PaintScale paintScale = new GrayPaintScale(minZ, maxZ);
        //PaintScale paintScale = new ColorPaintScale(minZ, maxZ);
        //PaintScale paintScale = new ColorPaintScale(minZ, maxZ, Color.red);
        PaintScale paintScale = new ColorPaintScale(minZ, maxZ, Color.red, true);

        renderer.setPaintScale(paintScale);

        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setAxisOffset(new RectangleInsets(5, 5, 5, 5));
        //plot.setDomainGridlinesVisible(false);

        final JFreeChart chart = new JFreeChart(title, null, plot, false);
        chart.setBackgroundPaint(Color.white);
        //chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 0, 1000, Color.green));

        PaintScaleLegend psl = new PaintScaleLegend(paintScale, zAxis);
        psl.setAxisOffset(5.0);
        //psl.setPosition(RectangleEdge.BOTTOM);
        psl.setPosition(RectangleEdge.RIGHT);
        psl.setMargin(new RectangleInsets(10, 10, 10, 10));
        chart.addSubtitle(psl);

        return(chart);

    }

    private JFreeChart createContourChart(String title, String xAxisLabel, String yAxisLabel, String zAxisLabel, boolean legend, ContourDataset dataset) {

        //final java.lang.String title = "Contour Plot";
        //final java.lang.String xAxisLabel = "X Values";
        //final java.lang.String yAxisLabel = "Y Values";
        //final java.lang.String zAxisLabel = "Color Values";

        ValueAxis xAxis = null;    /** The x-axis. */
        NumberAxis yAxis = null;   /** The y-axis. */
        ColorBar zColorBar = null; /** The z-axis. */
        //boolean zIsVertical = false;    /** Flag for vertical z-axis. */
        boolean xIsDate = false;  /** Flag for x is date axis. */
        boolean xIsLog = false;    /** Flag for x is log. */
        boolean yIsLog = false;    /** Flag for y is log. */
        boolean zIsLog = false;    /** Flag for z is log. */
        boolean xIsInverted = false;    /** Flag for x is inverted. */
        boolean yIsInverted = false;    /** Flag for y is inverted. */
        boolean zIsInverted = false;    /** Flag for z is inverted. */
        double ratio = 0.0;     /** The ratio. */

        if (xIsDate) {
            xAxis = new DateAxis(xAxisLabel);
            xIsLog = false; // force axis to be linear when displaying dates
        }
        else {
            if (xIsLog) {
                xAxis = new LogarithmicAxis(xAxisLabel);
            }
            else {
                xAxis = new NumberAxis(xAxisLabel);
            }
        }

        if (yIsLog) {
            yAxis = new LogarithmicAxis(yAxisLabel);
        }
        else {
            yAxis = new NumberAxis(yAxisLabel);
        }

        if (zIsLog) {
            zColorBar = new ColorBar(zAxisLabel);
        }
        else {
            zColorBar = new ColorBar(zAxisLabel);
        }

        if (xAxis instanceof NumberAxis) {
            ((NumberAxis) xAxis).setAutoRangeIncludesZero(false);
            ((NumberAxis) xAxis).setInverted(xIsInverted);
        }

        yAxis.setAutoRangeIncludesZero(false);

        yAxis.setInverted(yIsInverted);

        if (!xIsDate) {
            ((NumberAxis) xAxis).setLowerMargin(0.0);
            ((NumberAxis) xAxis).setUpperMargin(0.0);
        }

        yAxis.setLowerMargin(0.0);
        yAxis.setUpperMargin(0.0);

        zColorBar.getAxis().setInverted(zIsInverted);
        zColorBar.getAxis().setTickMarksVisible(true);

        final ContourPlot plot = new ContourPlot(dataset, xAxis, yAxis, zColorBar);

        if (xIsDate) {
            ratio = Math.abs(ratio); // don't use plot units for ratios when x axis is date
        }
        plot.setDataAreaRatio(ratio);

        final JFreeChart chart = new JFreeChart(null, null, plot, false);
        chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 0, 1000, Color.green));

        return chart;

    }

    /**
     * The starting point for the demo.
     *
     * @param args  ignored.
     *
     * @throws IOException  if there is a problem saving the file.
     */
    public static void main(String[] args) throws IOException, Exception {
        BioJavaChart chart = new BioJavaChart();
        chart.commandLine(args);
        JFrame frame = new JFrame();
        frame.setContentPane(chart);
        frame.pack();
        RefineryUtilities.centerFrameOnScreen(frame);
        frame.setVisible(true);

        WindowListener windowListener = new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
                    System.exit(0);
                }
            };
        frame.addWindowListener(windowListener);

    }

    private String[] getParams(String[] argsList, int index) {
        int counter = 0;
        while ((index+counter < argsList.length) && (!argsList[index+counter].startsWith("-"))) {
            counter++;
        }
        String[] params = new String[counter];
        for(int i=0; i < counter; i++) {
            params[i] = argsList[index+i];
        }
        //out.println("params_end="+params[params.length-1]);
        return(params);
    }


    public void commandLine(String[] args) throws IOException, Exception{
        String chartType = null;
        String fileName = null;
        int bins = -1;
        String title = null;
        String xAxisLabel = null;
        String yAxisLabel = null;
        String zAxisLabel = null;
        boolean legend = true;
        String[] fileNames = null;
        String saveFile = null;
        int saveWidth = 600;
        int saveHeigth = 400;
        String transposeString = "No";
        boolean transpose = false;

        if (args.length < 2 || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-?")) {
            out.println("\nBioJavaChart - The BioJava Charting Tool (using jFreeChart)\n");

            out.println(" Usage: BioJavaChart  <Options>\n");
            out.println(" Chart Options:");
            out.println();
            out.println("\t -hist \t\t <Values-File> <Values-File> ...");
            out.println();
            out.println("\t\t\t (one value per line)");
            out.println();
            out.println("\t -lines \t <X tab Y1 tab Y2... Values File>");
            out.println();
            out.println("\t\t\t (two or more values per line)");
            out.println();
            out.println("\t -scatter \t <X tab Y1 tab Y2... Values File>");
            out.println();
            out.println("\t\t\t (two or more values per line)");
            out.println();
            out.println("\t -heatmap \t <X tab Y tab Z Values File>");
            out.println();
            out.println("\t\t\t (three values per line)");
            out.println();
            out.println("\t -surface \t <X tab Y tab Z Values File>");
            out.println();
            out.println("\t\t\t (three values per line)");
            out.println();
            out.println("\t -bar \t\t <label tab value File>");
            out.println("\t\t\t (two values per line)");
            out.println();
            out.println("\t\t\t <label tab mean tab stdev File>");
            out.println("\t\t\t (three values per line)");
            out.println();
            out.println("\t -roc \t\t <fileName tab DataName File>");
            out.println();
            out.println("\t\t\t (two values per line)");
            out.println();
            out.println(" Other Options:");
            out.println();
            out.println("\t -bins \t\t <Num-Bins>");
            out.println();
            out.println("\t\t\t (if histogram, then default bins = MAX)");
            out.println("\t\t\t (if ROC curve, then default bins = "+maxBins+")");
            out.println();
            //out.println("\n\t -rocStepSize \t <Threshold-Step-Size> (default .005)");
            //out.println();
            out.println("\t -title \t <Chart-Title>");
            out.println();
            out.println("\t -xLabel \t <x-Label>");
            out.println();
            out.println("\t -yLabel \t <y-Label>");
            out.println();
            out.println("\t -zLabel \t <z-Label>");
            out.println();
            out.println("\t -legend \t <Yes/No> (default Yes)");
            out.println();
            //out.println("\t -transpose \t <Yes/No> (default "+transposeString+")");
            out.println("\t -transpose \t Use the transpose of the data matrix");
            out.println();
            out.println("\t -save \t\t <fileName> <width> <heigth>");
            out.println("\n\t\t\t (save format is png)");
            out.println("\t\t\t (example: -save hello.png 600 400)");
            out.println();
            out.println(" Data Files:");
            out.println();
            out.println("\t DatasetLabels \t Are set by the first commented line in");
            out.println("\t\t\t the values file:");
            out.println();
            out.println("\t example:\t # Y1-Label tab Y2-Label tab Y3-Label...");
            out.println();
            System.exit(-1);
        }

        for (int a=0; a < args.length; a++) {
            if  (args[a].equalsIgnoreCase("-hist")) {
                chartType = new String("histogram");
                fileNames = getParams(args,a+1);
            }
            else if  (args[a].equalsIgnoreCase("-bins")) {
                bins = Integer.parseInt(args[a+1]);
                this.maxBins = bins;
            }
            else if  (args[a].equalsIgnoreCase("-lines")) {
                chartType = new String("lines");
                fileName = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-scatter")) {
                chartType = new String("scatter");
                fileName = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-heatmap")) {
                chartType = new String("heatmap");
                fileName = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-surface")) {
                chartType = new String("surface");
                fileName = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-bar")) {
                chartType = new String("bar");
                fileNames = getParams(args,a+1);
                //fileName = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-roc")) {
                chartType = new String("roc");
                fileName = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-title")) {
                //title = getParams(args,a+1);
                title = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-xLabel")) {
                //xAxisLabel = getParams(args,a+1);
                xAxisLabel = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-yLabel")) {
                //yAxisLabel = getParams(args,a+1);
                yAxisLabel = new String(args[a+1]);
            }
            else if  (args[a].equalsIgnoreCase("-zLabel")) {
                //zAxisLabel = getParams(args,a+1);
                zAxisLabel = new String(args[a+1]);
            }
            else if (args[a].equals("-transpose")) {
                // transposeString = new String(args[a+1]);
                // transpose = StringTools.parseBoolean(transposeString);
                transpose = true;
            }
            else if  (args[a].equalsIgnoreCase("-save")) {
                //zAxisLabel = getParams(args,a+1);
                saveFile = new String(args[a+1]);
                saveWidth = Integer.parseInt(args[a+2]);;
                saveHeigth = Integer.parseInt(args[a+3]);;
            }
            else if  (args[a].equalsIgnoreCase("-legend")) {
                if (args[a+1].equalsIgnoreCase("Yes")) {
                    legend = true;
                }
                else {
                    legend = false;
                }
            }
        }

        //String path = getPathFromPathname(fileName);
        String path = null;

        if (chartType.equalsIgnoreCase("histogram")) {
            //String[] lines = FileTools.readStrings(fileName, 0, null);
            //String[] datasetLabels = getDatasetLabels(lines[0], fileName);
            double[][] values = getValuesMatrix(fileNames);
            initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, null, legend, values, bins);
        }
        else if (chartType.equalsIgnoreCase("lines")) {
            String[] lines = FileTools.readStrings(fileName, 0, null);
            //double[][] values = convertLinesToTwoDoubles(lines);
            String[] datasetLabels = getDatasetLabels(lines[0], fileName);
            double[][] values = StringTools.toDoubleMatrix(lines, "#", false);
            initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, values, bins);
        }
        else if (chartType.equalsIgnoreCase("scatter")) {
            String[] lines = FileTools.readStrings(fileName, 0, null);
            //double[][] values = convertLinesToTwoDoubles(lines);
            String[] datasetLabels = getDatasetLabels(lines[0], fileName);
            double[][] values = StringTools.toDoubleMatrix(lines, "#", false);
            initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, values, bins);
        }
        else if (chartType.equalsIgnoreCase("bar")) {
            CategoryDataset dataset = null;

            for (int i=0; i < fileNames.length; i++) {
                String[] lines = FileTools.readStrings(fileNames[i], 0, null);
                String[] datasetLabels = getDatasetLabels(lines[0], fileNames[i]);
                String[][] values = StringTools.toStringMatrix(lines, "#", false);
                // If just one value, then not whisker plot
                if (values.length ==2) {
                    if (dataset == null){
                        dataset = new DefaultCategoryDataset();
                    }
                    for (int j=0; j < values[0].length; j++) {
                        ((DefaultCategoryDataset)dataset).addValue(Double.parseDouble(values[1][j]), datasetLabels[0], values[0][j]);
                    }
                }
                // If 2 values, then we have a mean and stdev, use a whisker plot
                else {
                    if (dataset == null){
                        dataset = new DefaultStatisticalCategoryDataset();;
                    }
                    for (int j=0; j < values[0].length; j++) {
                        ((DefaultStatisticalCategoryDataset)dataset).add(Double.parseDouble(values[1][j]), Double.parseDouble(values[2][j]), datasetLabels[0], values[0][j]);
                    }
                }
            }
            chart = createBarChart(title, xAxisLabel, yAxisLabel, legend, dataset);
            createDefaultChartPanel(chart);
            postInitialize();

            //initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, legend, dataset);
        }

        else if (chartType.equalsIgnoreCase("heatmap") || chartType.equalsIgnoreCase("surface")) {
            String[] lines = FileTools.readStrings(fileName, 0, null);
            String[] datasetLabels = getDatasetLabels(lines[0], fileName);

            //double[][] valuesMatrix = convertLinesToThreeDoubles(lines);
            //datasetLabels = getDatasetLabels(lines[0], fileName);
            //double[][] valuesMatrix = StringTools.toDoubleMatrix(lines);

            // values are in the table format X1 Y1 Z1; X2 Y2 Z3; ....
            String[][] labelsMatrix = StringTools.toStringMatrix(lines, "#", false);

            //out.println("labelsMatrix.length="+labelsMatrix.length); == 3!!!!!!!!

            if (transpose) {
                String[][] newLabelsMatrix = new String[labelsMatrix.length][labelsMatrix[0].length];

                for (int i=0; i < labelsMatrix[0].length; i++) {
                    //swap X and Y
                    // newLabelsMatrix[i][0] = labelsMatrix[i][1];
                    // newLabelsMatrix[i][1] = labelsMatrix[i][0];
                    newLabelsMatrix[0][i] = labelsMatrix[1][i];
                    newLabelsMatrix[1][i] = labelsMatrix[0][i];
                    newLabelsMatrix[2][i] = labelsMatrix[2][i];
                }

                labelsMatrix = newLabelsMatrix;
            }

            double[][] valuesMatrix = new double[labelsMatrix.length][labelsMatrix[0].length];
            LinkedHashMap[] axisLabels = new LinkedHashMap[3];

            // create the LinkedHashMaps for any x, y, or z column that contains labels (not numbers)
            for (int j=0; j < 3; j++) {
                if (!MathTools.isNumber(labelsMatrix[j][0])) {
                    axisLabels[j] = new LinkedHashMap();
                }
                else {
                    axisLabels[j] = null;
                }

            }

            // Go though the matrix, filling in the valuesMatrix with numbers and maintaining sets of labels
            for (int i=0; i < labelsMatrix[0].length; i++) {
                for (int j=0; j < labelsMatrix.length; j++) {

                    if (axisLabels[j] != null) {
                        Double value = (Double) axisLabels[j].get(labelsMatrix[j][i]);

                        if (value == null) {
                            value = new Double(axisLabels[j].size());
                            axisLabels[j].put(labelsMatrix[j][i], value);
                        }

                        valuesMatrix[j][i] = value.doubleValue();
                    }
                    else { // the labels are numbers so just copy them
                        valuesMatrix[j][i] = Double.parseDouble(labelsMatrix[j][i]);
                    }
                }
            }

            // create string[]'s from the LinkedHashMaps of labels
            String[][] axisStrings = new String[3][];
            for (int j=0; j < 3; j++) {
                if (axisLabels[j] != null) {
                    axisStrings[j] = (String[]) axisLabels[j].keySet().toArray(new String[0]);
                }
                else {
                    axisStrings[j] = null;
                }
            }

            initialize(chartType, title, xAxisLabel, axisStrings[0], yAxisLabel, axisStrings[1], zAxisLabel, axisStrings[2], datasetLabels, legend, valuesMatrix, bins);
        }

        else if (chartType.equalsIgnoreCase("roc")) {
            String[] lines = FileTools.readStrings(fileName, 0, null);
            if (bins == -1) {
                bins = maxBins;
            }
            String[][] modelStrings = StringTools.toStringMatrix(lines, "#", false);
            fileNames = modelStrings[0];
            String[] datasetLabels = modelStrings[1];

            initialize(chartType, title, xAxisLabel, yAxisLabel, zAxisLabel, datasetLabels, legend, path, fileNames, bins);
        }

        // Save the chart to file if needed
        if (saveFile != null) {
            final ChartRenderingInfo info = new ChartRenderingInfo(new StandardEntityCollection());
            final File file1 = new File(saveFile);
            ChartUtilities.saveChartAsPNG(file1, chart, saveWidth, saveHeigth, info);
        }
    }


    private String format(double num, int minDecimalPlaces, int maxDecimalPlaces) {
        NumberFormat format = NumberFormat.getInstance();
        format.setMinimumFractionDigits(minDecimalPlaces);
        format.setMaximumFractionDigits(maxDecimalPlaces);
        return(format.format(num));
    }


    private double[][] getValuesMatrix(String[] fileNames) {
        double[][] values = new double[fileNames.length][];
        for (int i=0; i < values.length; i++) {
            String[] lines = FileTools.readStrings(fileNames[i], 0, null);
            values[i] = StringTools.toDoubles(lines, "#");
        }
        return(values);
    }

    private double[][][] getValuesMatrices(String[] fileNames) {
        double[][][] values = new double[fileNames.length][][];
        for (int i=0; i < values.length; i++) {
            String[] lines = FileTools.readStrings(fileNames[i], 0, null);
            values[i] = StringTools.toDoubleMatrix(lines, "#", true);
        }
        return(values);
    }

    // If there are no dataset labels in the first line of the datafile, then use the fileName as the dataset label
    public String[] getDatasetLabels(String lineString, String fileName) {
        if (lineString.trim().startsWith("#")) {
            //out.println(lineString);
            lineString = (lineString.trim().replaceFirst("#","")).trim();
            //out.println(lineString);
            return(lineString.split("\\t+"));
        }
        String[] stringArray = {fileName};
        return(stringArray);
    }

    public void createDefaultChartPanel(JFreeChart chart) {
        chartPanel = new ChartPanel(chart);
        setLayout(new BorderLayout());
        add(chartPanel, BorderLayout.CENTER);
    }

    public void createPanScrollZoomPanel(JFreeChart chart) {
        createPanScrollZoomPanel(chart, false);
    }

    public void createPanScrollZoomPanel(JFreeChart chart, boolean includeBinsSlider) {

        this.toolBar = createToolbar(includeBinsSlider);
        setLayout(new BorderLayout());
        add(this.toolBar, BorderLayout.SOUTH);

        //final JFreeChart chart = createChart();

        //this.scrollBar.setModel(new DefaultBoundedRangeModel());
        //recalcScrollBar(chart.getPlot());
        this.scrollBar.getModel().addChangeListener(this);

        this.chartPanel = new ChartPanel(chart) {
                public void autoRangeBoth() {
                    recalcScrollBar();
                    //out.println("Use 'Fit all' button");
                }
            };

        chart.addChangeListener(this);

        // enable zoom
        actionPerformed(new ActionEvent(this, 0, ACTION_CMD_ZOOM_BOX));

        // MouseListeners for pan function
        this.chartPanel.addMouseListener(this);
        this.chartPanel.addMouseMotionListener(this);

        // remove popup menu to allow panning
        // with right mouse pressed
        // this.chartPanel.setPopupMenu(null);

        //this.scrollBar.setModel(new DefaultBoundedRangeModel());
        recalcScrollBar(chart.getPlot());
        add(this.chartPanel);
    }

    /**
     * Creates a sample chart.
     *
     * @return a sample chart.
     */
    private JFreeChart createSampleChart() {

        final XYSeriesCollection primaryJFreeColl = new XYSeriesCollection();
        final XYSeries left1 = new XYSeries("Left 1");
        left1.add(1, 2);
        left1.add(2.8, 5.9);
        left1.add(3, null);
        left1.add(3.4, 2);
        left1.add(5, -1);
        left1.add(7, 1);
        primaryJFreeColl.addSeries(left1);

        final XYSeriesCollection secondaryJFreeColl = new XYSeriesCollection();
        final XYSeries right1 = new XYSeries("Right 1");
        right1.add(3.5, 2.2);
        right1.add(1.2, 1.3);
        right1.add(5.7, 4.1);
        right1.add(7.5, 7.4);
        secondaryJFreeColl.addSeries(right1);

        final NumberAxis xAxis = new NumberAxis("X");
        xAxis.setAutoRangeIncludesZero(false);
        xAxis.setAutoRangeStickyZero(false);

        final NumberAxis primaryYAxis = new NumberAxis("Y1");
        primaryYAxis.setAutoRangeIncludesZero(false);
        primaryYAxis.setAutoRangeStickyZero(false);

        // create plot
        final XYItemRenderer y1Renderer = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        y1Renderer.setSeriesPaint(0, Color.blue);
        y1Renderer.setToolTipGenerator(new StandardXYToolTipGenerator());
        final XYPlot xyPlot = new XYPlot(primaryJFreeColl, xAxis, primaryYAxis, y1Renderer);

        // 2nd y-axis

        final NumberAxis secondaryYAxis = new NumberAxis("Y2");
        secondaryYAxis.setAutoRangeIncludesZero(false);
        secondaryYAxis.setAutoRangeStickyZero(false);

        xyPlot.setRangeAxis(1, secondaryYAxis);
        xyPlot.setDataset(1, secondaryJFreeColl);

        xyPlot.mapDatasetToRangeAxis(1, 1);
        xyPlot.mapDatasetToDomainAxis(1, 1);

        final XYItemRenderer y2Renderer = new StandardXYItemRenderer(
            StandardXYItemRenderer.SHAPES_AND_LINES
                                                                     );
        y2Renderer.setToolTipGenerator(new StandardXYToolTipGenerator());
        xyPlot.setRenderer(1, y2Renderer);

        // set some fixed y-dataranges and remember them
        // because default chartPanel.autoRangeBoth()
        // would destroy them

        ValueAxis axis = xyPlot.getRangeAxis();
        this.primYMinMax[0] = -5;
        this.primYMinMax[1] = 15;
        axis.setLowerBound(this.primYMinMax[0]);
        axis.setUpperBound(this.primYMinMax[1]);

        axis = xyPlot.getRangeAxis(1);
        this.secondYMinMax[0] = -1;
        this.secondYMinMax[1] = 10;
        axis.setLowerBound(this.secondYMinMax[0]);
        axis.setUpperBound(this.secondYMinMax[1]);

        // Title + legend

        final String title = "To pan in zoom mode hold right mouse pressed";
        final JFreeChart ret = new JFreeChart(title, null, xyPlot, true);
        final TextTitle textTitle = new TextTitle(
            "(but you can only pan if the chart was zoomed before)"
                                                  );
        ret.addSubtitle(textTitle);
        return ret;
    }

    /**
     * Creates the toolbar.
     *
     * @return the toolbar.
     */
    private JToolBar createToolbar(boolean includeBinsSlider) {
        final JToolBar toolbar = new JToolBar();
        toolbar.setLayout(new BorderLayout());

        JPanel nested1 = new JPanel(new BorderLayout());
        JPanel nested1Left = new JPanel(new FlowLayout(FlowLayout.LEADING, 1, 0)); // FlowLayout
        nested1.add(nested1Left, BorderLayout.WEST);
        toolbar.add(nested1, BorderLayout.CENTER);


        final ButtonGroup groupedButtons = new ButtonGroup();

        // ACTION_CMD_PAN
        this.panButton = new JToggleButton();
        prepareButton(this.panButton, ACTION_CMD_PAN, "Pan", "Click & drag to pan up, down, left, or right");
        groupedButtons.add(this.panButton);
        nested1Left.add(this.panButton);

        // ACTION_CMD_ZOOM_BOX
        this.zoomButton = new JToggleButton();
        prepareButton(this.zoomButton, ACTION_CMD_ZOOM_BOX, "Zoom", "Click & drag to create zoom-in rectangle");
        groupedButtons.add(this.zoomButton);
        this.zoomButton.setSelected(true); // no other makes sense after startup
        nested1Left.add(this.zoomButton);

        // end of toggle-button group for select/pan/zoom-box
        //nested1Left.addSeparator();

        // ACTION_CMD_ZOOM_IN
        this.zoomInButton = new JButton();
        prepareButton(this.zoomInButton, ACTION_CMD_ZOOM_IN, "+", "Zoom in");
        nested1Left.add(this.zoomInButton);

        // ACTION_CMD_ZOOM_OUT
        this.zoomOutButton = new JButton();
        prepareButton(this.zoomOutButton, ACTION_CMD_ZOOM_OUT, "-", "Zoom out");
        nested1Left.add(this.zoomOutButton);

        // ACTION_CMD_ZOOM_TO_FIT
        this.fitButton = new JButton();
        prepareButton(this.fitButton, ACTION_CMD_ZOOM_TO_FIT, "Reset", "Reset pan and zoom");
        nested1Left.add(this.fitButton);

        //nested1Left.addSeparator();

        this.scrollBar = new JScrollBar(JScrollBar.HORIZONTAL);
        // scrollBar.setPreferredSize(new Dimension(0, 25));
        // int ht = (int) zoomButton.getPreferredSize().getHeight();
        // scrollBar.setMaximumSize(new Dimension(0, ht));
        // scrollBar.setPreferredSize(new Dimension(0, ht));
        this.scrollBar.setModel(new DefaultBoundedRangeModel());

        nested1.add(this.scrollBar, BorderLayout.CENTER);

        JPanel nested2 = null;
        if (includeBinsSlider) {
            nested2 = new JPanel(new BorderLayout());
            nested2.add(new JLabel("Number of bins:   1"), BorderLayout.WEST);
            nested2.add(new JLabel(""+maxBins), BorderLayout.EAST);
            binsSlider = new JSlider(SwingConstants.HORIZONTAL);
            nested2.add(binsSlider, BorderLayout.CENTER);

            binsSlider.setMinimum(1);
            binsSlider.setMaximum(maxBins);
            binsSlider.setValue(maxBins);

            binsSlider.addChangeListener(new SliderListener());

            SliderToolTips.enableSliderToolTips(binsSlider);

        }

        JPanel nested3 = null;
        if (includeLogLinear) {
            // nested3 = new JPanel(new FlowLayout(FlowLayout.CENTER,1,1));
            nested3 = new JPanel();
            nested3.setLayout(new BoxLayout(nested3, BoxLayout.X_AXIS));

            nested3.setPreferredSize(nested1.getPreferredSize());

            //this.linearXButton = new JRadioButton();
            this.linearXButton = new JToggleButton();
            linearXButton.setAlignmentX(Component.CENTER_ALIGNMENT);
            prepareButton(this.linearXButton, ACTION_CMD_X_LINEAR, "Linear Scale X-axis", "Use a linear scale for the x-axis");

            //this.logXButton = new JRadioButton();
            this.logXButton = new JToggleButton();
            logXButton.setAlignmentX(Component.CENTER_ALIGNMENT);
            prepareButton(this.logXButton, ACTION_CMD_X_LOG, "Log Scale X-axis", "Use a log scale for the x-axis");

            //this.linearYButton = new JRadioButton();
            this.linearYButton = new JToggleButton();
            linearYButton.setAlignmentX(Component.CENTER_ALIGNMENT);
            prepareButton(this.linearYButton, ACTION_CMD_Y_LINEAR, "Linear Scale Y-axis", "Use a linear scale for the y-axis");

            // this.logYButton = new JRadioButton();
            this.logYButton = new JToggleButton();
            logYButton.setAlignmentX(Component.CENTER_ALIGNMENT);
            prepareButton(this.logYButton, ACTION_CMD_Y_LOG, "Log Scale Y-axis", "Use a log scale for the y-axis");

            //this.unityLineCheckBox = new JCheckBox();
            this.unityLineCheckBox = new JToggleButton();
            unityLineCheckBox.setAlignmentX(Component.CENTER_ALIGNMENT);
            prepareButton(this.unityLineCheckBox, ACTION_CMD_UNITY_LINE, "Unity Line", "Draws a unity (y = x) line in the plot");
            unityLineCheckBox.setSelected(false);

            ButtonGroup xGroup = new ButtonGroup();
            xGroup.add(linearXButton);
            xGroup.add(logXButton);
            linearXButton.setSelected(true);

            ButtonGroup yGroup = new ButtonGroup();
            yGroup.add(linearYButton);
            yGroup.add(logYButton);
            linearYButton.setSelected(true);

            //////////////////////////////////////////////////////// - create layout
            nested3.add(Box.createHorizontalGlue());
            nested3.add(linearXButton);
            nested3.add(Box.createRigidArea(new Dimension(1,0)));
            nested3.add(logXButton);

            //nested3.addSeparator();
            nested3.add(Box.createRigidArea(new Dimension(20,0)));

            nested3.add(linearYButton);
            nested3.add(Box.createRigidArea(new Dimension(1,0)));
            nested3.add(logYButton);

            //nested3.addSeparator();
            nested3.add(Box.createRigidArea(new Dimension(20,0)));

            nested3.add(unityLineCheckBox);
            nested3.add(Box.createHorizontalGlue());
            //////////////////////////////////////////////////////// - end create layout
        }

        if (includeBinsSlider && !includeLogLinear) {
            toolbar.add(nested2, BorderLayout.SOUTH);
        }
        else if (!includeBinsSlider && includeLogLinear) {
            toolbar.add(nested3, BorderLayout.SOUTH);
        }
        else { // both!
            JPanel nested4 = new JPanel();
            nested4.setLayout(new BoxLayout(nested4, BoxLayout.Y_AXIS));
            nested4.add(nested2);
            nested4.add(nested3);
            toolbar.add(nested4, BorderLayout.SOUTH);
        }

        this.zoomOutButton.setEnabled(false);
        this.fitButton.setEnabled(false);
        this.scrollBar.setEnabled(false);

        toolbar.setFloatable(false);
        return toolbar;
    }

    public void drawUnityLine() {
        // XYLineAnnotation unityLineAnno = new XYLineAnnotation(0, 0, 1, 1, new BasicStroke(1f), Color.orange);

        Line2D.Double unityLine = new Line2D.Double(0, 0, 1, 1);
        XYShapeAnnotation unityLineAnno = new XYShapeAnnotation(unityLine, new BasicStroke(1f), Color.orange);

        chartPanel.getChart().getXYPlot().addAnnotation(unityLineAnno);

        unityLineCheckBox.setSelected(true);
    }

    /**
     * Prepares a button.
     *
     * @param button  the button.
     * @param actionKey  the action key.
     * @param buttonLabelText  the button label.
     * @param toolTipText  the tooltip text.
     */
    private void prepareButton(final AbstractButton button,
        final String actionKey,
        final String buttonLabelText,
        final String toolTipText) {
        // todo
        // as this action is empty and the button text is
        // redefined later, it can be safely removed ...
        //        Action action = new AbstractAction(actionKey) {
        //            public void actionPerformed(ActionEvent evt) {
        //                // ignored
        //            }
        //        };
        //        button.addActionListener(action);
        button.setActionCommand(actionKey);
        button.setText(buttonLabelText);
        button.setToolTipText(toolTipText);
        button.addActionListener(this);
    }

    /**
     * Sets the pan mode.
     *
     * @param val  a boolean.
     */
    private void setPanMode(final boolean val) {

        //this.chartPanel.setHorizontalZoom(!val);
        this.chartPanel.setDomainZoomable(!val);

        // chartPanel.setHorizontalAxisTrace(! val);

        //this.chartPanel.setVerticalZoom(!val);
        this.chartPanel.setRangeZoomable(!val);

        // chartPanel.setVerticalAxisTrace(! val);

        if (val) {
            //this.chartPanel.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            this.chartPanel.setCursor(Cursor.getPredefinedCursor(Cursor.MOVE_CURSOR));
        }
        else {
            //this.chartPanel.setCursor(Cursor.getDefaultCursor());
            this.chartPanel.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
        }
    }

    /**
     * Handles an action event.
     *
     * @param evt
     *            the event.
     */
    public void actionPerformed(final ActionEvent evt) {
        try {
            final String acmd = evt.getActionCommand();

            if (acmd.equals(ACTION_CMD_ZOOM_BOX)) {
                setPanMode(false);
            }
            else if (acmd.equals(ACTION_CMD_PAN)) {
                setPanMode(true);
            }
            else if (acmd.equals(ACTION_CMD_ZOOM_IN)) {
                final ChartRenderingInfo info = this.chartPanel.getChartRenderingInfo();
                final Rectangle2D rect = info.getPlotInfo().getDataArea();
                zoomBoth(rect.getCenterX(), rect.getCenterY(), ZOOM_FACTOR);
            }
            else if (acmd.equals(ACTION_CMD_ZOOM_OUT)) {
                final ChartRenderingInfo info = this.chartPanel.getChartRenderingInfo();
                final Rectangle2D rect = info.getPlotInfo().getDataArea();
                zoomBoth(rect.getCenterX(), rect.getCenterY(), 1 / ZOOM_FACTOR);
            }
            else if (acmd.equals(ACTION_CMD_ZOOM_TO_FIT)) {

                // X-axis (has no fixed borders)
                //this.chartPanel.autoRangeHorizontal();
                this.chartPanel.restoreAutoDomainBounds();

                //this.chartPanel.restoreAutoRangeBounds();

                // Y-Axes (autoRangeVertical)
                // not useful because of fixed borders
                final Plot plot = this.chartPanel.getChart().getPlot();

                ValueAxis rangeAxis = ((XYPlot)plot).getRangeAxis();
                rangeAxis.setAutoRange(true);

                // out.println("Lower Bound: "+rangeAxis.getLowerBound()+"\n");
                // out.println("Upper Bound: "+rangeAxis.getUpperBound()+"\n");
                this.primYMinMax[0] = rangeAxis.getLowerBound();
                this.primYMinMax[1] = rangeAxis.getUpperBound();

                rangeAxis.setAutoRange(false);
                numBinsChanged = false;

                if (plot instanceof ValueAxisPlot) {

                    final XYPlot vvPlot = (XYPlot) plot;
                    ValueAxis axis = vvPlot.getRangeAxis();
                    if (axis != null) {
                        axis.setLowerBound(this.primYMinMax[0]);
                        axis.setUpperBound(this.primYMinMax[1]);
                    }
                    if (plot instanceof XYPlot) {
                        final XYPlot xyPlot = (XYPlot) plot;
                        axis = xyPlot.getRangeAxis(1);
                        if (axis != null) {
                            axis.setLowerBound(this.secondYMinMax[0]);
                            axis.setUpperBound(this.secondYMinMax[1]);
                        }
                    }
                }
            }
            else if (acmd.equals(ACTION_CMD_X_LINEAR)) {
                final NumberAxis domainAxis = new NumberAxis("Relative Affinity (Ka)");
                chartPanel.getChart().getXYPlot().setDomainAxis(domainAxis);
            }
            else if (acmd.equals(ACTION_CMD_X_LOG)) {
                final NumberAxis domainAxis = new LogarithmicAxis("Relative Affinity (Ka)");
                // final ValueAxis domainAxis = new LogAxis("Relative Affinity (Ka)");
                chartPanel.getChart().getXYPlot().setDomainAxis(domainAxis);
                // drawUnityLine();
            }
            else if (acmd.equals(ACTION_CMD_Y_LINEAR)) {
                final NumberAxis rangeAxis = new NumberAxis("Relative Affinity (Ka)");
                chartPanel.getChart().getXYPlot().setRangeAxis(rangeAxis);
            }
            else if (acmd.equals(ACTION_CMD_Y_LOG)) {
                final NumberAxis rangeAxis = new LogarithmicAxis("Relative Affinity (Ka)");
                // final ValueAxis rangeAxis = new LogAxis("Relative Affinity (Ka)");
                chartPanel.getChart().getXYPlot().setRangeAxis(rangeAxis);
                // drawUnityLine();
            }
            else if (acmd.equals(ACTION_CMD_UNITY_LINE)) {
                if (unityLineCheckBox.isSelected()) {
                    drawUnityLine();
                }
                else {
                    chartPanel.getChart().getXYPlot().clearAnnotations();
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Handles a {@link ChangeEvent} (in this case, coming from the scrollbar).
     *
     * @param event  the event.
     */
    public void stateChanged(final ChangeEvent event) {
        try {
            final Object src = event.getSource();
            final BoundedRangeModel scrollBarModel = this.scrollBar.getModel();
            if (src == scrollBarModel) {
                final int val = scrollBarModel.getValue();
                final int ext = scrollBarModel.getExtent();

                final Plot plot = this.chartPanel.getChart().getPlot();
                if (plot instanceof XYPlot) {

                    //out.println("Hi!\n");

                    final XYPlot hvp = (XYPlot) plot;
                    final ValueAxis axis = hvp.getDomainAxis();

                    // avoid problems
                    this.chartPanel.getChart().removeChangeListener(this);

                    axis.setRange(val / this.scrollFactor, (val + ext) / this.scrollFactor);

                    // restore chart listener
                    this.chartPanel.getChart().addChangeListener(this);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Handles a {@link ChartChangeEvent}.
     *
     * @param event  the event.
     */
    public void chartChanged(final ChartChangeEvent event) {

        //out.println("Hi5!\n");

        try {

            if (event.getChart() == null) {
                return;
           }

            final BoundedRangeModel scrollBarModel = this.scrollBar.getModel();
            if ((scrollBarModel == null) && (this.binsSlider == null)) {
                //out.println("Hi3!\n");
                return;
            }

            boolean chartIsZoomed = false;

            final Plot plot = event.getChart().getPlot();
            if (plot instanceof XYPlot) {
                //if (true) {

                //out.println("Hi2!\n");

                final XYPlot hvp = (XYPlot) plot;
                final ValueAxis xAxis = hvp.getDomainAxis();
                final Range xAxisRange = xAxis.getRange();

                // avoid recursion
                scrollBarModel.removeChangeListener(this);

                // Added this to reset the scrollbar
                //recalcScrollBar(plot);

                final int low = (int) (xAxisRange.getLowerBound() * this.scrollFactor);
                scrollBarModel.setValue(low);
                final int ext = (int) (xAxisRange.getUpperBound() * this.scrollFactor - low);
                scrollBarModel.setExtent(ext);

                // restore
                scrollBarModel.addChangeListener(this);

                // check if zoomed horizontally
                //Range hdr = hvp.getHorizontalDataRange(xAxis);
                final Range hdr = hvp.getDataRange(xAxis);

                final double len = hdr == null ? 0 : hdr.getLength();
                chartIsZoomed |= xAxisRange.getLength() < len;
            }

            if (!chartIsZoomed && plot instanceof XYPlot) {
                // check if zoomed vertically
                final XYPlot vvp = (XYPlot) plot;
                ValueAxis yAxis = vvp.getRangeAxis();
                if (yAxis != null) {
                    chartIsZoomed = yAxis.getLowerBound() > this.primYMinMax[0]
                        || yAxis.getUpperBound() < this.primYMinMax[1];

                    // right y-axis
                    if (!chartIsZoomed && plot instanceof XYPlot) {
                        final XYPlot xyPlot = (XYPlot) plot;
                        yAxis = xyPlot.getRangeAxis(1);
                        if (yAxis != null) {
                            chartIsZoomed = yAxis.getLowerBound() > this.secondYMinMax[0]
                                || yAxis.getUpperBound() < this.secondYMinMax[1];
                        }
                    }
                }
            }

            // enable "zoom-out-buttons" if chart is zoomed
            // otherwise disable them
            this.panButton.setEnabled(chartIsZoomed);
            this.zoomOutButton.setEnabled(chartIsZoomed);
            this.scrollBar.setEnabled(chartIsZoomed);

            if (this.binsSlider == null) {
                this.fitButton.setEnabled(chartIsZoomed);
            }
            else {
                this.fitButton.setEnabled(chartIsZoomed || numBinsChanged);
            }

            if (!chartIsZoomed) {
                setPanMode(false);
                this.zoomButton.setSelected(true);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    // Mouse[Motion]Listeners for pan

    /**
     * Handles a mouse pressed event (to start panning).
     *
     * @param event  the event.
     */
    public void mousePressed(final MouseEvent event) {
        try {
            //if (this.panButton.isSelected() || this.panButton.isEnabled() && SwingUtilities.isRightMouseButton(event))
            if (this.panButton.isSelected() && !SwingUtilities.isRightMouseButton(event))
            {
                //final Rectangle2D dataArea = this.chartPanel.getScaledDataArea();
                final Rectangle2D dataArea = this.chartPanel.getScreenDataArea();
                final Point2D point = event.getPoint();
                if (dataArea.contains(point)) {
                    setPanMode(true);
                    this.panStartPoint = point;
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Handles a mouse released event (stops panning).
     *
     * @param event  the event.
     */
    public void mouseReleased(final MouseEvent event) {
        try {
            this.panStartPoint = null; // stop panning
            if (!this.panButton.isSelected()) {
                setPanMode(false);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Handles a mouse dragged event to perform panning.
     *
     * @param event  the event.
     */
    public void mouseDragged(final MouseEvent event) {
        try {
            if (this.panStartPoint != null) {
                //final Rectangle2D scaledDataArea = this.chartPanel.getScaledDataArea();
                final Rectangle2D scaledDataArea = this.chartPanel.getScreenDataArea();

                //this.panStartPoint = RefineryUtilities.getPointInRectangle(
                this.panStartPoint = getPointInRectangle(
                    this.panStartPoint.getX(),
                    this.panStartPoint.getY(),
                    scaledDataArea
                                                         );
                //final Point2D panEndPoint = RefineryUtilities.getPointInRectangle(
                final Point2D panEndPoint = getPointInRectangle(
                    event.getX(), event.getY(), scaledDataArea
                                                                );

                // horizontal pan

                final Plot plot = this.chartPanel.getChart().getPlot();
                if (plot instanceof XYPlot) {
                    final XYPlot hvp = (XYPlot) plot;
                    final ValueAxis xAxis = hvp.getDomainAxis();

                    if (xAxis != null) {
                        final double translatedStartPoint = xAxis.java2DToValue(
                            (float) this.panStartPoint.getX(),
                            scaledDataArea,
                            hvp.getDomainAxisEdge()
                                                                                );
                        final double translatedEndPoint = xAxis.java2DToValue(
                            (float) panEndPoint.getX(),
                            scaledDataArea,
                            hvp.getDomainAxisEdge()
                                                                              );
                        final double dX = translatedStartPoint - translatedEndPoint;

                        final double oldMin = xAxis.getLowerBound();
                        final double newMin = oldMin + dX;

                        final double oldMax = xAxis.getUpperBound();
                        final double newMax = oldMax + dX;

                        // do not pan out of range
                        if (newMin >= hvp.getDataRange(xAxis).getLowerBound()
                            && newMax <= hvp.getDataRange(xAxis).getUpperBound()) {
                            xAxis.setLowerBound(newMin);
                            xAxis.setUpperBound(newMax);
                        }
                    }
                }

                // vertical pan (1. Y-Axis)

                if (plot instanceof XYPlot) {
                    final XYPlot vvp = (XYPlot) plot;
                    final ValueAxis yAxis = vvp.getRangeAxis();

                    if (yAxis != null) {
                        final double translatedStartPoint = yAxis.java2DToValue(
                            (float) this.panStartPoint.getY(),
                            scaledDataArea,
                            vvp.getRangeAxisEdge()
                                                                                );
                        final double translatedEndPoint = yAxis.java2DToValue(
                            (float) panEndPoint.getY(),
                            scaledDataArea,
                            vvp.getRangeAxisEdge()
                                                                              );
                        final double dY = translatedStartPoint - translatedEndPoint;

                        final double oldMin = yAxis.getLowerBound();
                        final double newMin = oldMin + dY;

                        final double oldMax = yAxis.getUpperBound();
                        final double newMax = oldMax + dY;

                        // do not pan out of range
                        if (newMin >= this.primYMinMax[0] && newMax <= this.primYMinMax[1]) {
                            yAxis.setLowerBound(newMin);
                            yAxis.setUpperBound(newMax);
                        }
                    }
                }

                // vertical pan (2. Y-Axis)

                if (plot instanceof XYPlot) {
                    final XYPlot xyPlot = (XYPlot) plot;
                    final ValueAxis yAxis = xyPlot.getRangeAxis(1);

                    if (yAxis != null) {
                        final double translatedStartPoint = yAxis.java2DToValue(
                            (float) this.panStartPoint.getY(),
                            scaledDataArea,
                            xyPlot.getRangeAxisEdge(1)
                                                                                );
                        final double translatedEndPoint = yAxis.java2DToValue(
                            (float) panEndPoint.getY(),
                            scaledDataArea,
                            xyPlot.getRangeAxisEdge(1)
                                                                              );
                        final double dY = translatedStartPoint - translatedEndPoint;

                        final double oldMin = yAxis.getLowerBound();
                        final double newMin = oldMin + dY;

                        final double oldMax = yAxis.getUpperBound();
                        final double newMax = oldMax + dY;

                        if (newMin >= this.secondYMinMax[0] && newMax <= this.secondYMinMax[1]) {
                            yAxis.setLowerBound(newMin);
                            yAxis.setUpperBound(newMax);
                        }
                    }
                }

                // for the next time
                this.panStartPoint = panEndPoint;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Handles a mouse clicked event, in this case by ignoring it.
     *
     * @param event  the event.
     */
    public void mouseClicked(final MouseEvent event) {
        // ignored
    }

    /**
     * Handles a mouse moved event, in this case by ignoring it.
     *
     * @param event  the event.
     */
    public void mouseMoved(final MouseEvent event) {
        // ignored
    }

    /**
     * Handles a mouse entered event, in this case by ignoring it.
     *
     * @param event  the event.
     */
    public void mouseEntered(final MouseEvent event) {
        // ignored
    }

    /**
     * Handles a mouse exited event, in this case by ignoring it.
     *
     * @param event  the event.
     */
    public void mouseExited(final MouseEvent event) {
        // ignored
    }

    /**
     * Starting point for the demo.
     *
     * @param args  the command line arguments (ignored).
     */
    //     public static void main(final String[] args) {

    //         try {
    //             final String lookAndFeelClassName = WindowsLookAndFeel.class.getName();
    //             UIManager.setLookAndFeel(lookAndFeelClassName);
    //         }
    //         catch (Exception ex) {
    //             out.println(ex.getMessage());
    //         }

    //         final PanScrollZoomDemo demo = new PanScrollZoomDemo("Pan & Scroll & Zoom - Demo");
    //         demo.pack();
    //         demo.setVisible(true);

    //     }

    // PRIVATE

    public void recalcScrollBar() {
        recalcScrollBar(this.chartPanel.getChart().getPlot());
    }


    /**
     * Recalculates the scrollbar settings.
     *
     * @param plot  the plot.
     */
    public void recalcScrollBar(final Plot plot) {

        if (plot instanceof ContourPlot) {
            //out.println("Hi\n");
            return;
        }

        ValueAxis rangeAxis = ((XYPlot)plot).getRangeAxis();

        //this.chartPanel.autoRangeHorizontal();
        //this.chartPanel.autoRangeVertical();
        this.chartPanel.restoreAutoDomainBounds();
        this.chartPanel.restoreAutoRangeBounds();

        //out.println("Lower Bound: "+rangeAxis.getLowerBound()+"\n");
        //out.println("Upper Bound: "+rangeAxis.getUpperBound()+"\n");
        this.primYMinMax[0] = rangeAxis.getLowerBound();
        this.primYMinMax[1] = rangeAxis.getUpperBound();
        //this.primYMinMax[0] = 0;
        //this.primYMinMax[1] = 900;
        //rangeAxis.setLowerBound(this.primYMinMax[0]);
        //rangeAxis.setUpperBound(this.primYMinMax[1]);


        if (plot instanceof XYPlot) {
            final XYPlot hvp = (XYPlot) plot;
            final ValueAxis axis = hvp.getDomainAxis();

            axis.setLowerMargin(0);
            axis.setUpperMargin(0);

            final Range rng = axis.getRange();

            final BoundedRangeModel scrollBarModel = this.scrollBar.getModel();

            //added removeChangeListener
            scrollBarModel.removeChangeListener(this);

            final int len = scrollBarModel.getMaximum() - scrollBarModel.getMinimum();
            if (rng.getLength() > 0) {
                this.scrollFactor = len / rng.getLength();
            }

            final double dblow = rng.getLowerBound();
            final int ilow = (int) (dblow * this.scrollFactor);
            scrollBarModel.setMinimum(ilow);
            final int val = ilow;
            scrollBarModel.setValue(val);

            final double dbup = rng.getUpperBound();
            final int iup = (int) (dbup * this.scrollFactor);
            scrollBarModel.setMaximum(iup);
            final int ext = iup - ilow;
            scrollBarModel.setExtent(ext);

            scrollBarModel.addChangeListener(this);

            // Fit

            // X-axis (has no fixed borders)
            //this.chartPanel.autoRangeHorizontal();

            // Y-Axes) (autoRangeVertical
            // not useful because of fixed borders
            if (plot instanceof ValueAxisPlot) {
                //if (true) {

                //out.println("Hi4!\n");

                final XYPlot vvPlot = (XYPlot) plot;
                //ValueAxis rangeAxis = vvPlot.getRangeAxis();
                if (rangeAxis != null) {
                    //out.println("Hi6!\n");
                    rangeAxis.setLowerBound(this.primYMinMax[0]);
                    rangeAxis.setUpperBound(this.primYMinMax[1]);
                }
                if (plot instanceof XYPlot) {
                    final XYPlot xyPlot = (XYPlot) plot;
                    rangeAxis = xyPlot.getRangeAxis(1);
                    //out.println("Hi7!\n");
                    if (rangeAxis != null) {
                        //out.println("Hi8!\n");
                        rangeAxis.setLowerBound(this.secondYMinMax[0]);
                        rangeAxis.setUpperBound(this.secondYMinMax[1]);
                    }
                }
            }


        }
    }

    /**
     * Zooms in on an anchor point (measured in Java2D coordinates).
     *
     * @param x  the x value.
     * @param y  the y value.
     * @param zoomFactor  the zoomFactor < 1 == zoom in; else out.
     */
    private void zoomBoth(final double x, final double y, final double zoomFactor) {
        zoomHorizontal(x, zoomFactor);
        zoomVertical(y, zoomFactor);
    }

    /**
     * Decreases the range on the horizontal axis, centered about a Java2D x coordinate.
     * <P>
     * The range on the x axis is multiplied by zoomFactor
     *
     * @param x  the x coordinate in Java2D space.
     * @param zoomFactor  the zoomFactor < 1 == zoom in; else out.
     */
    private void zoomHorizontal(final double x, final double zoomFactor) {

        final JFreeChart chart = this.chartPanel.getChart();
        final ChartRenderingInfo info = this.chartPanel.getChartRenderingInfo();
        if (chart.getPlot() instanceof XYPlot) {
            final XYPlot hvp = (XYPlot) chart.getPlot();
            final ValueAxis axis = hvp.getDomainAxis();
            if (axis != null) {
                final double anchorValue = axis.java2DToValue(
                    (float) x, info.getPlotInfo().getDataArea(), hvp.getDomainAxisEdge()
                                                              );
                if (zoomFactor < 1.0) {
                    axis.resizeRange(zoomFactor, anchorValue);
                }
                else if (zoomFactor > 1.0) {
                    final Range range = hvp.getDataRange(axis);
                    adjustRange(axis, range, zoomFactor, anchorValue);
                }
            }
        }
    }

    /**
     * Decreases the range on the vertical axis, centered about a Java2D y coordinate.
     * <P>
     * The range on the y axis is multiplied by zoomFactor
     *
     * @param y  the y coordinate in Java2D space.
     * @param zoomFactor  the zoomFactor < 1 == zoom in; else out.
     */
    private void zoomVertical(final double y, final double zoomFactor) {

        final JFreeChart chart = this.chartPanel.getChart();
        final ChartRenderingInfo info = this.chartPanel.getChartRenderingInfo();

        // 1. (left) Y-Axis

        if (chart.getPlot() instanceof XYPlot) {
            final XYPlot vvp = (XYPlot) chart.getPlot();
            final ValueAxis primYAxis = vvp.getRangeAxis();
            if (primYAxis != null) {
                final double anchorValue =
                    primYAxis.java2DToValue(
                        (float) y, info.getPlotInfo().getDataArea(), vvp.getRangeAxisEdge()
                                            );
                if (zoomFactor < 1.0) {
                    // zoom in
                    primYAxis.resizeRange(zoomFactor, anchorValue);

                }
                else if (zoomFactor > 1.0) {
                    // zoom out
                    final Range range = new Range(this.primYMinMax[0], this.primYMinMax[1]);
                    adjustRange(primYAxis, range, zoomFactor, anchorValue);
                }
            }

            // 2. (right) Y-Axis

            if (chart.getPlot() instanceof XYPlot) {
                final XYPlot xyp = (XYPlot) chart.getPlot();
                final ValueAxis secYAxis = xyp.getRangeAxis(1);
                if (secYAxis != null) {
                    final double anchorValue =
                        secYAxis.java2DToValue(
                            (float) y,
                            info.getPlotInfo().getDataArea(),
                            xyp.getRangeAxisEdge(1));
                    if (zoomFactor < 1.0) {
                        // zoom in
                        secYAxis.resizeRange(zoomFactor, anchorValue);

                    }
                    else if (zoomFactor > 1.0) {
                        // zoom out
                        final Range range = new Range(this.secondYMinMax[0], this.secondYMinMax[1]);
                        adjustRange(secYAxis, range, zoomFactor, anchorValue);
                    }
                }
            }
        }
    }

    /**
     * used for zooming
     *
     * @param axis  the axis.
     * @param range  the range.
     * @param zoomFactor  the zoom factor.
     * @param anchorValue  the anchor value.
     */
    private void adjustRange(final ValueAxis axis, final Range range, final double zoomFactor,
        final double anchorValue) {

        if (axis == null || range == null) {
            return;
        }

        final double rangeMinVal = range.getLowerBound()
            - range.getLength() * axis.getLowerMargin();
        final double rangeMaxVal = range.getUpperBound()
            + range.getLength() * axis.getUpperMargin();
        final double halfLength = axis.getRange().getLength() * zoomFactor / 2;
        double zoomedMinVal = anchorValue - halfLength;
        double zoomedMaxVal = anchorValue + halfLength;
        double adjMinVal = zoomedMinVal;
        if (zoomedMinVal < rangeMinVal) {
            adjMinVal = rangeMinVal;
            zoomedMaxVal += rangeMinVal - zoomedMinVal;
        }
        double adjMaxVal = zoomedMaxVal;
        if (zoomedMaxVal > rangeMaxVal) {
            adjMaxVal = rangeMaxVal;
            zoomedMinVal -= zoomedMaxVal - rangeMaxVal;
            adjMinVal = Math.max(zoomedMinVal, rangeMinVal);
        }

        final Range adjusted = new Range(adjMinVal, adjMaxVal);
        axis.setRange(adjusted);
    }

    /**
     * Returns a point based on (x, y) but constrained to be within the bounds of a given
     * rectangle.
     *
     * @param x  the x-coordinate.
     * @param y  the y-coordinate.
     * @param area  the constraining rectangle.
     *
     * @return a point within the rectangle.
     */
    public static Point2D getPointInRectangle(double x, double y, final Rectangle2D area) {

        x = Math.max(area.getMinX(), Math.min(x, area.getMaxX()));
        y = Math.max(area.getMinY(), Math.min(y, area.getMaxY()));
        return new Point2D.Double(x, y);

    }

    private class SliderListener implements javax.swing.event.ChangeListener
    {
        public void stateChanged(javax.swing.event.ChangeEvent ce)
        {
            JSlider source = (JSlider) ce.getSource();
            if (! source.getValueIsAdjusting())
            {

                //out.println("Hi!\n");

                // val is between 1 and 100
                numBinsChanged = true;

                int binsSliderValue = source.getValue();

                // Have to re-create the dataset from scratch?!!
                //IntervalXYDataset dataset = createHistDataset(histDatasetLabels, histValues, binsSliderValue);

                if (histDataset != null) {
                    histDataset = createHistDataset(histDatasetLabels, histValues, binsSliderValue);

                    XYPlot xyplot = (XYPlot)(chart.getPlot());
                    xyplot.setDataset(histDataset);

                    ((HistogramDataset)histDataset).setType(HistogramType.RELATIVE_FREQUENCY);
                }

                else if (rocDataset != null) {
                    rocDataset.removeAllSeries();

                    RocData[] rocDataArray = new RocData[rocFileNames.length];
                    for (int i=0; i < rocFileNames.length; i++) {
                        if (rocPath != null) {
                            rocDataArray[i] = new RocData(rocDatasetLabels[i], rocPath + File.separator + rocFileNames[i], binsSliderValue);
                        }
                        else {
                            rocDataArray[i] = new RocData(rocDatasetLabels[i], rocFileNames[i], binsSliderValue);
                        }
                        String aucString = ": A.U.C. = "+format(rocDataArray[i].auc, 3, 3);

                        XYSeries series = null;
                        if (rocDatasetLabels != null) {
                            series = createXYSeries(rocDatasetLabels[i]+aucString, rocDataArray[i].rocArray[0], rocDataArray[i].rocArray[1]);
                        }
                        else {
                            series = createXYSeries(aucString, rocDataArray[i].rocArray[0], rocDataArray[i].rocArray[1]);
                        }

                        rocDataset.addSeries(series);
                    }
                }

                // chartPanel.repaint();


                // // X-axis (has no fixed borders)
                // //this.chartPanel.autoRangeHorizontal();
                // chartPanel.restoreAutoDomainBounds();

                // // Y-Axes (autoRangeVertical)
                // // not useful because of fixed borders
                // final Plot plot = chartPanel.getChart().getPlot();
                // if (plot instanceof ValueAxisPlot) {

                //     final XYPlot vvPlot = (XYPlot) plot;
                //     ValueAxis axis = vvPlot.getRangeAxis();
                //     if (axis != null) {
                //         axis.setLowerBound(primYMinMax[0]);
                //         axis.setUpperBound(primYMinMax[1]);
                //     }
                //     if (plot instanceof XYPlot) {
                //         final XYPlot xyPlot = (XYPlot) plot;
                //         axis = xyPlot.getRangeAxis(1);
                //         if (axis != null) {
                //             axis.setLowerBound(secondYMinMax[0]);
                //             axis.setUpperBound(secondYMinMax[1]);
                //         }
                //     }
                // }


            }
        }
    }

}
