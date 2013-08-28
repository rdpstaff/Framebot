/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

import java.io.*;

/**
 *
 * @author wangqion
 */
public class FrameBotStatIterator {
    private FrameBotReaderCore core;
    
    public FrameBotStatIterator(File infile, boolean ignoreAlignment) throws FileNotFoundException, IOException{
        InputStreamReader is = new InputStreamReader(new FileInputStream(infile));
        char firstChar = (char) is.read();
        if ( firstChar == '>'){
            is.close();
            core = new FrameBotAlignOutputReader(infile, ignoreAlignment);
        }else {
            is.close();
            core = new FrameBotStatLineReader(infile);
        }
    }
    
    public boolean hasNext(){
        return core.hasNext();
    }
    
    public FrameBotStat next(){
        return core.next();
    }
    
    
    public void close(){
        core.close();
    }
    
}
