# Coding Style

- Avoid runtime OOP at all costs. Do not introduce inheritance-based polymorphism,
  virtual dispatch, owning base-class pointers, or framework-style object hierarchies.
- Compile-time polymorphism is acceptable. Prefer templates, concepts, overloads,
  value types, and plain data structures when behavior needs to vary by type.
- Prefer pure free functions over member functions.
- If a member function exists, keep it simple: getters, setters, trivial convenience
  wrappers, or direct delegation to an accompanying pure free function.
- Any nontrivial operation should be available as a pure free function even when a
  member-function convenience wrapper is also present.
