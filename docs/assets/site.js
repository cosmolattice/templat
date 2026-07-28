/* ═══════════════════════════════════════════════════════════════════════
   TempLat website — shared behaviour: scroll reveal, hero logo draw-on,
   sticky nav (+ hero parallax) and code-block copy buttons.

   Every piece is optional: each looks for its own hooks and does nothing
   when they are absent, so the same file serves a page with a hero
   (index.html) and one without (vocabulary.html).

   Companion stylesheet: site.css
   ═══════════════════════════════════════════════════════════════════════ */
(function () {
  'use strict';
  var reduced = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  /* ── scroll reveal ─────────────────────────────────────────────────
     Content must never end up permanently invisible, so this is built
     defensively: IntersectionObserver when available, a scroll-position
     fallback when it is not, and a timer that force-reveals everything
     if neither has fired. Worst case the animation is skipped — the
     page is still readable.                                          */
  var items = Array.prototype.slice.call(document.querySelectorAll('.reveal'));

  function showAll() { items.forEach(function (el) { el.classList.add('in'); }); }

  if (reduced) {
    showAll();
  } else {
    /* stagger siblings inside a grid so cards cascade rather than snap in together */
    var seen = new Map();
    items.forEach(function (el) {
      var p = el.parentNode, n = seen.get(p) || 0;
      seen.set(p, n + 1);
      el.style.setProperty('--d', (n % 6) * 0.065 + 's');
    });

    if ('IntersectionObserver' in window) {
      /* the observer delivers an initial callback for every target, visible
         or not — so this flag tells us it is alive, independently of whether
         anything has actually scrolled into view yet */
      var alive = false;
      var io = new IntersectionObserver(function (entries) {
        alive = true;
        entries.forEach(function (e) {
          if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); }
        });
      }, { rootMargin: '0px 0px -12% 0px', threshold: 0.08 });
      items.forEach(function (el) { io.observe(el); });

      /* failsafe: observer never delivered — show everything rather than nothing */
      setTimeout(function () { if (!alive) showAll(); }, 2500);
    } else {
      /* no observer: reveal on scroll position instead */
      var check = function () {
        var h = window.innerHeight;
        items = items.filter(function (el) {
          if (el.getBoundingClientRect().top < h * 0.88) { el.classList.add('in'); return false; }
          return true;
        });
      };
      window.addEventListener('scroll', check, { passive: true });
      window.addEventListener('resize', check);
      check();
    }
  }

  /* ── hero logo draw-on ─────────────────────────────── */
  var logo = document.getElementById('heroLogo');
  if (logo) requestAnimationFrame(function () { logo.classList.add('play'); });

  /* ── sticky nav + hero parallax ─────────────────────────────────────
     The nav is hidden by CSS only while JS is running (.js .nav), so it
     is this script's job to bring it back. On a page with a hero it
     slides in past the fold; on a page without one there is nothing to
     scroll past, so it is shown straight away and stays put.         */
  var nav = document.getElementById('nav');
  var hero = document.querySelector('.hero');
  var grid = document.getElementById('heroGrid');
  var ticking = false;

  if (nav && !hero) nav.classList.add('on');

  function onScroll() {
    var y = window.pageYOffset || document.documentElement.scrollTop;
    if (nav && hero) nav.classList.toggle('on', y > window.innerHeight * 0.72);
    if (grid && !reduced) grid.style.setProperty('--par', (y * 0.16) + 'px');
    ticking = false;
  }
  window.addEventListener('scroll', function () {
    if (!ticking) { ticking = true; requestAnimationFrame(onScroll); }
  }, { passive: true });
  onScroll();

  /* ── copy buttons ──────────────────────────────────── */
  document.querySelectorAll('.copy').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var pre = btn.closest('.code').querySelector('pre');
      var text = pre.innerText;

      function ok() {
        var was = btn.textContent;
        btn.textContent = 'Copied';
        btn.classList.add('done');
        setTimeout(function () { btn.textContent = was; btn.classList.remove('done'); }, 1600);
      }

      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(ok, legacy);
      } else {
        legacy();
      }

      /* file:// and older browsers have no async clipboard */
      function legacy() {
        var ta = document.createElement('textarea');
        ta.value = text;
        ta.setAttribute('readonly', '');
        ta.style.cssText = 'position:absolute;left:-9999px;top:0';
        document.body.appendChild(ta);
        ta.select();
        try { document.execCommand('copy'); ok(); } catch (e) { /* nothing else to try */ }
        document.body.removeChild(ta);
      }
    });
  });
})();
