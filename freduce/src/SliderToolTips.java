import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.awt.EventQueue;

import javax.swing.Action;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.JSlider;
import javax.swing.JComponent;
import java.text.NumberFormat;

// @author Santhosh Kumar T - santhosh@in.fiorano.com
// @author Todd Riley
public class SliderToolTips{

    public static void enableSliderToolTips(final JSlider slider) {
        enableSliderToolTips(slider, 1.0, 0, 0, 0);
    }

    public static void enableSliderToolTips(final JSlider slider, final double scalar, final double offset, final int minDecimalPlaces, final int maxDecimalPlaces) {

        slider.addChangeListener(new ChangeListener(){
                private boolean adjusting = false;
                private String oldTooltip;

                public void stateChanged(ChangeEvent e){

                    if (slider.getModel().getValueIsAdjusting()) {

                        // System.out.println("is adjusting");

                        if (!adjusting){
                            oldTooltip = slider.getToolTipText();
                            adjusting = true;
                        }

                        String toolTipText;
                        if (! ((scalar == 1.0) && (offset == 0)) ) {
                            toolTipText = format((slider.getValue() * scalar) + offset, minDecimalPlaces, maxDecimalPlaces);
                        }
                        else {
                            toolTipText = String.valueOf(slider.getValue());
                        }

                        slider.setToolTipText(toolTipText);
                        hideToolTip(slider); // to avoid flickering :)
                        postToolTip(slider);
                    }
                    else{
                        // System.out.println("is NOT adjusting");

                        hideToolTip(slider);
                        slider.setToolTipText(oldTooltip);
                        adjusting = false;
                        oldTooltip = null;
                    }
                }
            });
    }

    public static void enableSliderToolTips(final JSlider slider, final boolean normalize, final boolean reverse, final int minDecimalPlaces, final int maxDecimalPlaces) {

        slider.addChangeListener(new ChangeListener(){
                private boolean adjusting = false;
                private String oldTooltip;

                public void stateChanged(ChangeEvent e){

                    if (slider.getModel().getValueIsAdjusting()) {

                        // System.out.println("is adjusting");

                        if (!adjusting){
                            oldTooltip = slider.getToolTipText();
                            adjusting = true;
                        }

                        String toolTipText;
                        if (normalize && !reverse) {
                            toolTipText = format(
                                ((double)slider.getValue()) / slider.getMaximum(),
                                minDecimalPlaces,
                                maxDecimalPlaces);
                        }
                        else if (!normalize && reverse) {
                            toolTipText = String.valueOf(slider.getMaximum() - (slider.getValue() - slider.getMinimum()) );
                        }
                        else if (normalize && reverse) {
                            toolTipText = format(
                                (1.0) - ((((double)slider.getValue()) -  slider.getMinimum())/ slider.getMaximum()),
                                minDecimalPlaces,
                                maxDecimalPlaces);
                        }
                        else {
                            // why did you call this method?
                            toolTipText = String.valueOf(slider.getMaximum());
                        }

                        slider.setToolTipText(toolTipText);
                        hideToolTip(slider); // to avoid flickering :)
                        postToolTip(slider);
                    }
                    else{
                        // System.out.println("is NOT adjusting");

                        hideToolTip(slider);
                        slider.setToolTipText(oldTooltip);
                        adjusting = false;
                        oldTooltip = null;
                    }
                }
            });
    }

    private static String format(double num, int minDecimalPlaces, int maxDecimalPlaces) {
        NumberFormat format = NumberFormat.getInstance();
        format.setMinimumFractionDigits(minDecimalPlaces);
        format.setMaximumFractionDigits(maxDecimalPlaces);
        return(format.format(num));
    }


    private static boolean isInteger(double aNumber) {
        if (aNumber == Math.round(aNumber)) {
        //if (aNumber == Math.floor(aNumber)) {
        //if ((aNumber == Math.rint(aNumber)) && (!Double.isInfinite(aNumber))) {
            return true;
        }
        return false;
    }

    /*-------------------------------------------------[ Manual ToolTips ]---------------------------------------------------*/

    public static void postToolTip(JComponent comp) {
        Action action = comp.getActionMap().get("postTip");
        if (action == null) // no tooltip
            return;
        ActionEvent ae = new ActionEvent(comp, ActionEvent.ACTION_PERFORMED,
            "postTip", EventQueue.getMostRecentEventTime(), 0);
        action.actionPerformed(ae);
    }

    public static void hideToolTip(JComponent comp) {
        Action action = comp.getActionMap().get("hideTip");
        if (action == null) {// no tooltip
            return;
        }
        ActionEvent ae = new ActionEvent(comp, ActionEvent.ACTION_PERFORMED,
            "hideTip", EventQueue.getMostRecentEventTime(), 0);
        action.actionPerformed(ae);
    }
}
