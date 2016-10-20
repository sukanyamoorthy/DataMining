/*@author : sukanya moorthy
 *@version :1.0
 * 
 * Class represents a Gene .
 * Each gene is represented by number of sample_id's it belongs
 * to and list of sample_id's
 * 
 */

package org.dm.apriori;

import java.util.ArrayList;

public class Apriori {
	int count;
	ArrayList<String> listofsamples;

	public int getCount() {
		return count;
	}

	public void setCount(int count) {
		this.count = count;
	}

	public ArrayList<String> getListofsamples() {
		return listofsamples;
	}

	public void setListofsamples(ArrayList<String> listofsamples) {
		this.listofsamples = listofsamples;
	}

}
