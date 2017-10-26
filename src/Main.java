import org.apache.commons.math3.special.Erf;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import java.awt.*;
import java.util.HashMap;


public class Main extends ApplicationFrame {

    /**
     * Creates a new demo.
     *
     * @param title the frame title.
     */
    static HashMap<Double,Double> Mt = new HashMap<>();
    static HashMap<Double,Double> Et = new HashMap<>();
    public Main(final String title) {

        super(title);

        final XYDataset dataset = createDataset();
        final JFreeChart chart = createChart(dataset);
        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        setContentPane(chartPanel);

    }

    /**
     * Creates a sample dataset.
     *
     * @return a sample dataset.
     */
    private XYDataset createDataset() {

        final XYSeries series1 = new XYSeries("First");
        for (HashMap.Entry<Double,Double> doubleDoubleEntry : Mt.entrySet()) {
            series1.add(doubleDoubleEntry.getKey(),doubleDoubleEntry.getValue());
        }
        final XYSeries series2 = new XYSeries("Second");
        for (HashMap.Entry<Double,Double> doubleDoubleEntry : Et.entrySet()) {
            series2.add(doubleDoubleEntry.getKey(),doubleDoubleEntry.getValue());
        }

        final XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series1);
        dataset.addSeries(series2);

        return dataset;

    }

    /**
     * Creates a chart.
     *
     * @param dataset the data for the chart.
     * @return a chart.
     */
    private JFreeChart createChart(final XYDataset dataset) {

        // create the chart...
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "Model",      // chart title
                "Mt",                      // x axis label
                "c",                      // y axis label
                dataset,                  // data
                PlotOrientation.VERTICAL,
                true,                     // include legend
                true,                     // tooltips
                false                     // urls
        );

        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
        chart.setBackgroundPaint(Color.white);

//        final StandardLegend legend = (StandardLegend) chart.getLegend();
        //      legend.setDisplaySeriesShapes(true);

        // get a reference to the plot for further customisation...
        final XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.lightGray);
        //    plot.setAxisOffset(new Spacer(Spacer.ABSOLUTE, 5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);

        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(1, false);
        plot.setRenderer(renderer);

        // change the auto tick unit selection to integer units only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        // OPTIONAL CUSTOMISATION COMPLETED.

        return chart;

    }

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

    /**
     * Starting point for the demonstration application.
     *
     * @param args ignored.
     */
    public static void main(final String[] args) {
        tipomain();
        final Main demo = new Main("Model");
        demo.pack();
        RefineryUtilities.centerFrameOnScreen(demo);
        demo.setVisible(true);
        for (HashMap.Entry<Double,Double> doubleDoubleEntry : Et.entrySet()) {
            System.out.println(doubleDoubleEntry.getValue());
        }

    }



    static void tipomain() {
        System.out.println("Begin");
        double alpha = 1.0, betta = 1.0, M, D2, u = 10.0, c = 0.01, v = 0.0, t = 100.0;
        M = alpha / betta;
        D2 = (2.0 * alpha) / (Math.pow(betta, 2.0));
        //double D = Math.sqrt(D2);
        double x1 = 1.0;
        // ArrayList<Double> Mt = new ArrayList <Double>();


        while (c < 2.15) {
            double x2 = ((c * (t - v)) / (u + c * v)) + 1;

            double tmp1 = normal(((Math.sqrt(u + c * v)) / (c * Math.sqrt(D2) * Math.sqrt(x2))) * (x2 * (1 - c * M) - 1));
            double tmp3 = normal(((Math.sqrt(u + c * v)) / (c * Math.sqrt(D2) * Math.sqrt(x1))) * (x1 * (1 - c * M) - 1));
            double tmp2 = Math.exp(((2 * (u + c * v)) / (c * c * D2)) * (1 - c * M)) * normal(-(Math.sqrt(u + c * v)) / (c * Math.sqrt(D2 * x2)) * (x2 * (1 - c * M) + 1));
            double tmp4 = Math.exp(((2 * (u + c * v)) / (c * c * D2)) * (1 - c * M)) * normal(-(Math.sqrt(u + c * v)) / (c * Math.sqrt(D2 * x1)) * (x1 * (1 - c * M) + 1));
            double tmp5 = normal((-(Math.sqrt(u + c * v)) / (c * Math.sqrt(D2 * x2))) * (x2 * (1 - c * M) + 1));
            double tmp6 = normal((-(Math.sqrt(u + c * v)) / (c * Math.sqrt(D2 * x1))) * (x1 * (1 - c * M) + 1));
            double tmp = (tmp1 + tmp2) - (tmp3 + tmp4);

            double Cf = betta/(4*alpha*c);
            double Cs = betta/(4*alpha*c);
           /* double Ft = (-(c*c*D2)/(u+c*v))*((tmp1+tmp2)-(tmp3+tmp4)) + 2*(1-c*M)*Math.exp(((2 * (u + c * v)) / (c * c * D2)) * (1 - c * M))*
                    (tmp5-tmp6)-((2*c*Math.sqrt(D2))/Math.sqrt(Math.PI*x2*(u+c*v)))*Math.exp(-((u+c*v)/(2*x2*c*c*D2))*Math.pow(x2*(1-c*M)-1 ,2.0))+
                    ((2*c*Math.sqrt(D2))/Math.sqrt(Math.PI*x1*(u+c*v)))*Math.exp(-((u+c*v)/(2*x1*c*c*D2))*Math.pow(x1*(1-c*M)-1 ,2.0));*/
            double Ft = (-(c*c*D2)/(u+c*v))*((tmp1 + tmp2) - (tmp3 + tmp4)) + 2*(1 - c*M)*(tmp2 - tmp4)-
                    (((2*c*Math.sqrt(D2))/Math.sqrt(Math.PI*x1*(u+c*v)))*((Math.exp(((- (u + c*v)) / (2*x2*c*c*D2))*Math.pow(x2*(1 - c*M)-1 , 2.0))
                            -(Math.exp(((- (u + c*v)) / (2*x1*c*c*D2))*Math.pow(x1*(1 - c*M)-1 , 2.0))))));

            /*double St = (-(3*c*c*D2)/(u+c*v))*((tmp1+tmp2)-(tmp3+tmp4)) + 2*(1-c*M)*(3-4*((u+c*v)/(c*c*D2))*(1-c*M))*Math.exp(((2 * (u + c * v)) / (c * c * D2)) * (1 - c * M))*
                    (tmp5-tmp6)-((Math.sqrt(2.0)*c*Math.sqrt(D2))/(Math.sqrt(Math.PI*(u+c*v))*Math.pow(x2, 1.5))*(3*(1-((u+c*v)/(c*c*D2))*(1-c*M))*x2+(u+c*v)/(c*c*D2))*Math.exp((-((u+c*v)/(2*x2*c*c*D2)))*Math.pow(x2*(1-c*M)-1 ,2.0)))+
                    ((Math.sqrt(2.0)*c*Math.sqrt(D2))/(Math.sqrt(Math.PI*(u+c*v))*Math.pow(x1, 1.5))*(3*(1-((u+c*v)/(c*c*D2))*(1-c*M))*x1+(u+c*v)/(c*c*D2))*Math.exp((-((u+c*v)/(2*x2*c*c*D2)))*Math.pow(x1*(1-c*M)-1 ,2.0)));*/
            double St = (-(3*c*c*D2)/(u+c*v))* ((tmp1+tmp2)-(tmp3+tmp4)) + 2*(1-c*M)*(3-4*((u+c*v)/(c*c*D2))*(1-c*M))*(tmp2 - tmp4)-
                    (((Math.sqrt(2.0)*c*Math.sqrt(D2))/(Math.sqrt(Math.PI*(u+c*v))*Math.pow(x2, 1.5))*(3*(1-((u+c*v)/(c*c*D2))*(1-c*M))*x2+(u+c*v)/(c*c*D2)))*((Math.exp(((- (u + c*v)) / (2*x2*c*c*D2))*Math.pow(x2*(1 - c*M)-1 , 2.0))))
                            -((Math.sqrt(2.0)*c*Math.sqrt(D2))/(Math.sqrt(Math.PI*(u+c*v))*Math.pow(x1, 1.5))*(3*(1-((u+c*v)/(c*c*D2))*(1-c*M))*x1+(u+c*v)/(c*c*D2)))*((Math.exp(((- (u + c*v)) / (2*x1*c*c*D2))*Math.pow(x1*(1 - c*M)-1 , 2.0)))));

            Et.put(c,(tmp + Cf*Ft + Cs*St));
            Mt.put(c, tmp);
            c += 0.05;

        }



        System.out.println("end");

    }


    static double normal(double x) {

        //return (1.0/Math.sqrt(2.0*Math.PI)*integrator.integrate(1000,f, -Float.MAX_VALUE,x));
        return 0.5 * (1.0 + Erf.erf(x));
    }

}
