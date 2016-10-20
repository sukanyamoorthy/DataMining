/*@author : sukanya moorthy
 *@version :1.0
 * 
 * Util contains all the required functions used in calculating Apriori.
 *
 */

package org.dm.apriori;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Util {
	public static HashMap<String, Apriori> read_file(String filename) {
		// Parses file and generates a map with Gene and its related sample_id's
		HashMap<String, Apriori> mp = new HashMap<String, Apriori>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String key = "";
			try {
				String line = in.readLine();
				while (line != null) {
					String[] tokens = line.split("\t");
					int i = 1;
					for (String val : tokens) {
						if (val.equalsIgnoreCase("UP") || val.equalsIgnoreCase("Down")
								|| !val.toLowerCase().contains("Sample".toLowerCase())) {
							if (val.equalsIgnoreCase("UP") || (val.equalsIgnoreCase("Down"))) {
								key = val + i;
							} else
								key = val;
							if (mp.get(key) != null) {
								Apriori ob = mp.get(key);
								ob.setCount(ob.getCount() + 1);
								ArrayList<String> existingsamples = ob.getListofsamples();
								existingsamples.add(tokens[0]);
								ob.setListofsamples(existingsamples);
								mp.put(key, ob);
							} else {
								ArrayList<String> samples = new ArrayList<String>();
								samples.add(tokens[0]);
								Apriori ob = new Apriori();
								ob.setCount(1);
								ob.setListofsamples(samples);
								mp.put(key, ob);
							}
							i++;
						}
					}
					line = in.readLine();
				}
			} catch (IOException e) {

				e.printStackTrace();

			}
		} catch (FileNotFoundException e) {

			e.printStackTrace();
		}
		return mp;
	}

	public static void runApriori() {
		// Executes Apriori for the static file and rules
		HashMap<String, Apriori> file_map = Util.read_file("gene_expression.txt");
		HashMap<Set<String>, Integer> l1 = Util.generateL1(file_map);
		HashMap<Set<String>, Integer> prune = l1;
		int tot = 100;
		int support = 50;
		int support_val = (tot * support) / 100;
		prune = prune(prune, support_val);
		HashMap<String, Integer> rules = Util.associations(prune, file_map, support_val);
		System.out.println("########## Template 1 ##########");
		ArrayList<String> check = new ArrayList<String>();
		check.add("UP6");
		Util.query_evaluation(rules, check, "1", "Rule", true);
		check.add("UP1");
		check.add("Down10");
		Util.query_evaluation(rules, check, "1", "Rule", true);
	}

	public static HashMap<Set<String>, Integer> generateL1(HashMap<String, Apriori> freq_1) {
		// Generates candidate set
		HashMap<Set<String>, Integer> set_map = new HashMap<Set<String>, Integer>();
		for (String key : freq_1.keySet()) {
			Set<String> set = new HashSet<String>();
			set.add(key);
			set_map.put(set, freq_1.get(key).count);
		}
		return set_map;
	}

	public static HashMap<String, Integer> associations(HashMap<Set<String>, Integer> freq_1,
			HashMap<String, Apriori> mp, int val) {
		// Builds association based on confidence : 60
		int k = 2;
		HashMap<String, Integer> rules = new HashMap<String, Integer>();

		freq_1 = selfJoin(freq_1, k, mp);
		freq_1 = prune(freq_1, val);
		rules = rule_generation(freq_1, mp);
		freq_1 = selfJoin(freq_1, k + 1, mp);
		freq_1 = prune(freq_1, val);
		System.out.println("Total No of Associations " + rules.size());
		return rules;
	}

	public static void print_frequent_itemset(HashMap<String, Apriori> mp)
	// helper function to debug
	{
		int tot = 0;
		for (String key : mp.keySet()) {
			Apriori ob = mp.get(key);
			if (ob.count >= 50) {
				tot = tot + 1;
				// System.out.println(key+" \t"+ob.count+"\t");
			}

		}
		System.out.println("Count is: " + tot);

	}

	public static ArrayList<String> intersection(ArrayList<String> list1, ArrayList<String> list2) {
		// Returns intersecting elements between two lists
		ArrayList<String> list = new ArrayList<String>();

		for (String t : list1) {
			if (list2.contains(t)) {
				list.add(t);
			}
		}

		return list;
	}

	public static int intersection(Set<String> set, HashMap<String, Apriori> mp) {
		// Returns no of intersecting elements
		List list = new ArrayList(set);
		ArrayList<String> sample = new ArrayList<String>();
		sample = mp.get(list.get(0)).listofsamples;
		for (int i = 1; i < set.size(); i++) {
			// sample.retainAll(mp.get(list.get(i)).listofsamples);
			sample = intersection(sample, mp.get(list.get(i)).listofsamples);
		}
		return sample.size();
	}

	public static HashMap<Set<String>, Integer> selfJoin(HashMap<Set<String>, Integer> mp, int k,
			HashMap<String, Apriori> freq_1) {
		// Joins two sets
		HashMap<Set<String>, Integer> set_map = new HashMap<Set<String>, Integer>();
		List<Set<String>> key_set = new ArrayList<Set<String>>();
		key_set.addAll(mp.keySet());

		int len = key_set.size();
		for (int i = 0; i < len; i++) {
			Set<String> s1 = key_set.get(i);
			for (int j = i + 1; j < len; j++) {
				Set<String> s2 = key_set.get(j);
				Set<String> newset = new HashSet<String>();
				newset.addAll(s1);
				newset.addAll(s2);
				int count = 0;
				for (String x : s1) {
					if (s1.equals(s2)) {
						count++;
					}
				}
				if (newset.size() == k) {
					count = intersection(newset, freq_1);
					set_map.put(newset, count);
				}

			}
		}
		return set_map;
	}

	public static HashMap<Set<String>, Integer> prune(HashMap<Set<String>, Integer> mp, int val) {
		// Prunes the set based on support
		HashMap<Set<String>, Integer> result = new HashMap<Set<String>, Integer>();
		for (Set<String> s : mp.keySet()) {
			if (mp.get(s) >= val) {
				result.put(s, mp.get(s));
			}
		}
		return result;
	}

	public static HashMap<String, Integer> rule_generation(HashMap<Set<String>, Integer> freq_1,
			HashMap<String, Apriori> mp) {
		// Generates association rules based on confidence=60
		HashMap<String, Integer> set_map = new HashMap<String, Integer>();
		double con_n = 0;
		double con_d = 0;
		double con = 0;
		List<Set<String>> key_set = new ArrayList<Set<String>>();
		key_set.addAll(freq_1.keySet());
		for (int i = 0; i < key_set.size(); i++) {
			Set<String> s11 = key_set.get(i);
			con = freq_1.get(s11);
			List<String> s1 = new ArrayList<String>(s11);
			String tmp = "";
			String tmp2 = "";
			for (int j = 0; j < s1.size(); j++) {
				if (tmp.equals(""))
					tmp = tmp + s1.get(j);
				else
					tmp = tmp + "->" + s1.get(j);
			}
			String[] rev = tmp.split("->");
			con_n = mp.get(rev[0]).count;
			con_d = mp.get(rev[1]).count;
			tmp2 = rev[1] + "->" + rev[0];
			if ((con / con_n) * 100 >= 60) {
				int val = (int) ((con / con_n) * 100);
				set_map.put(tmp, val);
			}
			if ((con / con_d) * 100 >= 60) {
				int val = (int) ((con / con_d) * 100);
				set_map.put(tmp2, val);
			}

		}

		return set_map;
	}

	public static HashMap<String, Integer> query_evaluation(HashMap<String, Integer> set_map, ArrayList<String> rules,
			String type, String category, Boolean flag) {
		// Evaluates Query for the template "RULE HAS 1 OF ANY [INPUT LIST]"
		HashMap<String, Integer> result = new HashMap<String, Integer>();
		int cn = 0;
		if (category.equals("Rule")) {
			if (type.equals("1")) {
				for (String key : set_map.keySet()) {
					String[] token = key.split("->");
					List<String> head = Arrays.asList(token[0].split(","));
					List<String> body = Arrays.asList(token[1].split(","));

					Set<String> newSet = new HashSet<String>(head);

					newSet.addAll(body);
					head = new ArrayList<String>(newSet);
					int flg = 0;
					for (String t : head) {
						if (rules.contains(t)) {
							flg++;
							if (flg > 1) {
								break;
							}
						}
					}
					if (flg == 1) {
						result.put(key, 1);
						cn++;
					}
				}
			}

			if (flag)
			{
			System.out.println("\n"+category + " " + "HAS " + type + " OF " + rules.toString() + ": " + cn);
			System.out.println(result);
			}

		}

		return result;
	}

}
