/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

/**
 *
 * @author wangqion
 */
public abstract class FrameBotReaderCore {
    
    public abstract FrameBotStat next();
    public abstract boolean hasNext();
     
    public FrameBotStat getStatLine( String statline){
        String[] tokens = statline.trim().split("\\t");
        return new FrameBotStat(tokens[1], tokens[2], Integer.parseInt(tokens[3]), 
                Integer.parseInt(tokens[4]), Double.parseDouble(tokens[5]), 
                Integer.parseInt(tokens[6]), Integer.parseInt(tokens[7]), Boolean.parseBoolean(tokens[8]));
    }
    
    public FrameBotStat getAlignment( String[] statlines){
        String[] tokens = statlines[1].trim().split("\\t");
        return new FrameBotStat(tokens[1], tokens[2], Integer.parseInt(tokens[3]), 
                Integer.parseInt(tokens[4]), Double.parseDouble(tokens[5]), 
                Integer.parseInt(tokens[6]), Integer.parseInt(tokens[7]), Boolean.parseBoolean(tokens[8]), statlines);
    }
    
    
    public abstract void close();
    
}
