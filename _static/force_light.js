(function () {
  try {
    // Force light theme
    document.documentElement.dataset.theme = 'light';
    // Persist if theme toggle code looks at localStorage
    localStorage.setItem('theme', 'light');
  } catch (e) {}
})();
