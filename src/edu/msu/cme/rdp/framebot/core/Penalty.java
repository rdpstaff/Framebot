/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.core;

/**
 *
 * @author wangqion
 */
public class Penalty {

    public static final int gap_extension_cost = 1;
    public static final int v = 10;   // gap opening penalty
    public static final int u = -gap_extension_cost;   // extension cost
    public static final int FRAME_SHIFT_PENALTY = 10;
    public static final int w = getIndelCost(1);

    public static int getIndelCost(int k) {
        return -(gap_extension_cost * k + v);
    }
}
