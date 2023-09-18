""" 
Creator: Letizia D'Achille 

Functions:

	number2base_rep(n, b)
		Computes the representation of n in base b
	admissible(n, b)
		Determines if n is b-admissible
	count_admissible(b, start, end)
		Counts all b-admissible numbers within the range [start, end)
	count_admissible_width(b, width)
		Counts all b-admissible numbers whose b-representation has width digits
	largest_multi_admissible(L, start, end)
		Finds the largest number within the range [start, end) that is b-admissible for all b in L
"""

def number2base_rep(n, b):
	"""
	Computes the representation of n in base b

		Parameters:
			n (int): A decimal non-negative integer
			b (int): An integer 2 <= b <= 10
		
		Returns:
			s (string): String of the representation of n in base b
	"""
	s = str(n % b)
	n //= b
	while n > 0:
		s = str(n % b) + s
		n //= b
	return s

def admissible(n, b):
	"""
	Determines if n is b-admissible

		Parameters:
			n (int): A decimal non-negative integer
			b (int): An integer 2 <= b <= 10
		
		Returns:
			b (boolean): True if n is b-admissible, False otherwise
	"""
	s = number2base_rep(n, b)
	for i in range(1, len(s)//2 + 1):
		ind = i
		for j in range(0, len(s)-i):
			if s[j] == s[j+i]:
				if j == ind - 1:
					return False
			else:
				ind = j+i+1
	return True

def count_admissible(b, start, end):
	"""
	Counts all b-admissible numbers within the range [start, end)

		Parameters:
			b (int): An integer 2 <= b <= 10
			start (int): A non-negative integer
			end (int): A non-negative integer with start < end
		
		Returns:
			count (int): An integer with the number of b-admissible numbers within the range [start, end)
	"""
	count = 0
	for n in range(start, end):
		if admissible(n, b):
			count += 1
	return count

def count_admissible_width(b, width):
	"""
	Counts all b-admissible numbers whose b-representation has width digits

		Parameters:
			b (int): An integer 2 <= b <= 10
			width (int): A non-negative integer
		
		Returns:
			count (int): An integer with the number of b-admissible numbers whose b-representation has width digits
	"""
	return count_admissible(b, b**(width - 1) + (b**(width - 3) if width >= 3 else 0), b**width - (b**(width - 2) if width >= 2 else 0))

def largest_multi_admissible(L, start, end):
	"""
	Finds the largest number within the range [start, end) that is b-admissible for all b in L

		Parameters:
			L (list[int]): A list of integers 2 <= b <= 10
			start (int): A non-negative integer
			end (int): A non-negative integer with start < end
		
		Returns:
			n (int): An integer with the largest number within the range [start, end) that is b-admissible for all b in L. If no such number exists, returns None.
	"""
	for n in range(end - 1, start - 1, -1):
		flag = True
		for b in L:
			if not admissible(n, b):
				flag = False
				break
		if flag:
			return n
	return None