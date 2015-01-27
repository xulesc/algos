import std.stdio, std.algorithm, std.range, std.string;

dstring[] read_words(immutable dstring fileName) {
	dstring[] words;
	foreach(dchar[] w; lines(File("words.txt"))) {
		w = w.chomp().toLower();
		words ~= w.idup;
	}
	return words;	
}

dstring[] get_n_grams(immutable dstring word, immutable int size) {
	if(word.length < size)
		return [];
	int i;
	dstring[] ngrams;
	for(i = 0; i <= word.length - size; i++)
		ngrams ~= word[i..i+size];
	return ngrams;	
}

auto make_inverted_index(const dstring[] words) {
	dstring k = "||"; dstring[] v = ["||"];
	auto dict = [k:v];
	foreach(word; words) {
		foreach(ngram; get_n_grams(word, 1)) {
			dstring[] v1 = dict.get(ngram, []);
			v1 ~= word;
			dict[ngram] = v1;
		}
	}
	return dict;
}

void main() {
	// read the dictionary words
	dstring[] words = read_words("words.txt");
	// make inverted index
	auto dict = make_inverted_index(words);
}
