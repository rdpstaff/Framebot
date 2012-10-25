/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.Scanner;

/**
 *
 * @author wangqion
 */
public class FrameBotStatIterator implements Iterator{
    private Scanner scanner;
    
    public FrameBotStatIterator(File framebotResult) throws FileNotFoundException{
        scanner = new Scanner(framebotResult).useDelimiter(">");
    }
    
    public boolean hasNext(){
        return scanner.hasNext();
    }
    
    public FrameBotStat next(){
        String[] statline = scanner.next().split("\n");
        // the second line
        String[] tokens = statline[1].trim().split("\\t");
        return new FrameBotStat(tokens[1], tokens[2], Integer.parseInt(tokens[3]), 
                Integer.parseInt(tokens[4]),Integer.parseInt(tokens[5]), 
                Integer.parseInt(tokens[6]), Integer.parseInt(tokens[7]), Boolean.parseBoolean(tokens[8]));
    }
    
    public void remove(){
        throw new UnsupportedOperationException();
    }
    
    public void close(){
        scanner.close();
    }
    
}
