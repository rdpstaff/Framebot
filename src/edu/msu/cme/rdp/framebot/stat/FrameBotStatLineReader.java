/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author wangqion
 */
public class FrameBotStatLineReader extends FrameBotReaderCore{
    private BufferedReader reader = null;
    String statline = null;
    
    public FrameBotStatLineReader(File framebotResult) throws FileNotFoundException{
        reader = new BufferedReader(new FileReader(framebotResult));
    }
    
    public boolean hasNext(){
        try {
            statline = reader.readLine();
        } catch (IOException ex) {
            Logger.getLogger(FrameBotStatLineReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        return (statline != null);
    }
    
    public FrameBotStat next(){
        return getStatLine(statline);
    }
     
    public void close(){
        try {
            reader.close();
        } catch (IOException ex) {
            Logger.getLogger(FrameBotStatLineReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
